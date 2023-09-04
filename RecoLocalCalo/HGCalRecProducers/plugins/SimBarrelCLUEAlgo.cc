#include "RecoLocalCalo/HGCalRecProducers/plugins/SimBarrelCLUEAlgo.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"


#include "oneapi/tbb/task_arena.h"
#include "oneapi/tbb.h"
#include <limits>

using namespace hgcal_clustering;


template <typename T>
void SimBarrelCLUEAlgoT<T>::setThresholds(
                     edm::ESGetToken<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd> EBThres,
                     edm::ESGetToken<HcalPFCuts, HcalPFCutsRcd> HCALThres ){
  tok_ebThresholds_ = EBThres;
  tok_hcalThresholds_ = HCALThres;
}


template <typename T>
void SimBarrelCLUEAlgoT<T>::getEventSetupPerAlgorithm(const edm::EventSetup& es) {
  cells_.clear();
  numberOfClustersPerLayer_.clear();
  ebThresholds_ = &es.getData(tok_ebThresholds_);
  hcalThresholds_ = &es.getData(tok_hcalThresholds_);
  maxlayer_ = maxLayerIndex_; //for testing purposes
  cells_.resize(maxlayer_ + 1);
  numberOfClustersPerLayer_.resize(maxlayer_ + 1, 0);
}

template <typename T>
void SimBarrelCLUEAlgoT<T>::populate(const edm::PCaloHitContainer& hits) {
  std::unordered_map<uint32_t, double> hae;
  for (unsigned int i = 0; i < hits.size(); ++i) {
    PCaloHit hit = hits[i];
    std::unordered_map<uint32_t, double>::iterator it = hae.find(hit.id());
    if (it != hae.end()) {
      double energy = it->second + hit.energy();
      it->second = energy;
    } else {
      hae.emplace(hit.id(), hit.energy());
    }
  }
  //for (unsigned int i = 0; i < hae.size(); ++i) {
  for (auto hit = hae.begin(); hit != hae.end(); ++hit) {
    if (hit->second < 1e-2) continue;
    DetId detid(hit->first);
    if ((detid.det() == 3 || detid.det() == 4)) {
      if (detid.det() == 3 && detid.subdetId() <=5) {
	const GlobalPoint position = rhtools_.getPosition(detid);
	int layer = 0;

	if (detid.det() == DetId::Hcal) {
	  HcalDetId hid(detid);
	  layer = hid.depth();
	  if (detid.subdetId() == HcalSubdetector::HcalOuter) 
	    layer += 1;
	}
      
	cells_[layer].detid.push_back(detid);
	cells_[layer].eta.push_back(position.eta());
	cells_[layer].phi.push_back(position.phi());
	cells_[layer].weight.push_back(static_cast<float>(hit->second));
	cells_[layer].r.push_back(position.mag());
	float sigmaNoise = 0.f;
	if (detid.det() == DetId::Ecal) {
	  sigmaNoise = (*ebThresholds_)[detid];
	} else { 
	  const HcalPFCut* item = (*hcalThresholds_).getValues(detid);
	  sigmaNoise = item->seedThreshold();
	}
	cells_[layer].sigmaNoise.push_back(sigmaNoise);
      }
    }
  }
}

template <typename T>
void SimBarrelCLUEAlgoT<T>::prepareDataStructures(unsigned int l) {
  auto cellsSize = cells_[l].detid.size();
  cells_[l].rho.resize(cellsSize, 0.f);
  cells_[l].delta.resize(cellsSize, 9999999);
  cells_[l].nearestHigher.resize(cellsSize, -1);
  cells_[l].clusterIndex.resize(cellsSize);
  cells_[l].followers.resize(cellsSize);
  cells_[l].isSeed.resize(cellsSize, false);
  cells_[l].eta.resize(cellsSize, 0.f);
  cells_[l].phi.resize(cellsSize, 0.f);
  cells_[l].r.resize(cellsSize, 0.f);
}

template <typename T>
void SimBarrelCLUEAlgoT<T>::makeClusters() {
  tbb::this_task_arena::isolate([&] {
    tbb::parallel_for(size_t(0), size_t(maxlayer_ + 1), [&](size_t i) {
      prepareDataStructures(i);
      T lt;
      lt.clear();
      lt.fill(cells_[i].eta, cells_[i].phi);
      float delta_c;
      
      if constexpr (std::is_same_v<T, EBLayerTiles>)
	delta_c = vecDeltas_[0];
      else
	delta_c = vecDeltas_[1];
      float delta_r = vecDeltas_[2];
      
      calculateLocalDensity(lt, i, delta_c, delta_r);
      calculateDistanceToHigher(lt, i, delta_c, delta_r);
      numberOfClustersPerLayer_[i] = findAndAssignClusters(i, delta_c, delta_r);
      
      if (doSharing_)
	// Now running the sharing routine
	passSharedClusterIndex(lt, i, delta_c);
      });
    });
  for (unsigned int i = 0; i < maxlayer_ + 1; ++i) {
    setDensity(i);
  }
}

template <typename T>
std::vector<reco::BasicCluster> SimBarrelCLUEAlgoT<T>::getClusters(bool) {
  std::vector<int> offsets(numberOfClustersPerLayer_.size(), 0);

  int maxClustersOnLayer = numberOfClustersPerLayer_[0];

  if constexpr (!std::is_same_v<T, EBLayerTiles>) {
    for (unsigned layerId = 1; layerId < offsets.size(); ++layerId) {
      offsets[layerId] = offsets[layerId - 1] + numberOfClustersPerLayer_[layerId - 1];
      maxClustersOnLayer = std::max(maxClustersOnLayer, numberOfClustersPerLayer_[layerId]);
    }
  }

  auto totalNumberOfClusters = offsets.back() + numberOfClustersPerLayer_.back();
  clusters_v_.resize(totalNumberOfClusters);
  // Store the cellId and energy for each cluster
  // We need to store the energy as a pair here because we perform the energy splitting
  std::vector<std::vector<std::pair<int, float>>> cellsIdInCluster;
  cellsIdInCluster.reserve(maxClustersOnLayer);
  
  for (unsigned int layerId = 0; layerId < maxlayer_ + 1; ++layerId) {
    cellsIdInCluster.resize(numberOfClustersPerLayer_[layerId]);
    auto& cellsOnLayer = cells_[layerId];
    unsigned int numberOfCells = cellsOnLayer.detid.size();
    auto firstClusterIdx = offsets[layerId];

    for (unsigned int i = 0; i < numberOfCells; ++i) {
      //std::cout << "===> Cell " << i << "\n";
      auto clusterIndex = cellsOnLayer.clusterIndex[i];
      
      if (clusterIndex.size() == 1){
         cellsIdInCluster[clusterIndex[0]].push_back(
            std::make_pair(i, 1.f));
      }
      else if (clusterIndex.size() > 1) {
        std::vector<float> fractions (clusterIndex.size());
        
	for (unsigned int j = 0; j < clusterIndex.size(); j++) {
          const auto& seed = clusterIndex[j];
	  const auto& seedCell = cellsOnLayer.seedToCellIndex[seed];
	  const auto& seedEnergy = cellsOnLayer.weight[seedCell];
	  //std::cout << "Seed " << seed << std::endl;
	  //std::cout << "seedToCellIndex[seed]" << cellsOnLayer.seedToCellIndex[seed] << std::endl;
          // compute the distance
	  float dist = distance(i, seedCell, layerId) / T::type::cellWidthEta;
	  //std::cout << "Distance in cells unit " << dist << std::endl;
          fractions[j] = seedEnergy * std::exp(-(std::pow(dist,2))/(2*std::pow(T::type::showerSigma,2)));
	  //std::cout << "Cell " << i << " in cluster index " << j << " with fraction " << fractions[j] << std::endl;
        }
	//std::cout << "==================\n";
        auto tot_norm_fractions = std::accumulate(std::begin(fractions), std::end(fractions),  0.);

        for (unsigned int j = 0; j < clusterIndex.size(); j++) {
          float norm_fraction = fractions[j]/tot_norm_fractions;
	  //std::cout << "Adding Cell " << i << " to cluster " << clusterIndex[j] << " with energy " << cellsOnLayer.weight[i]*norm_fraction  << " (fraction " << norm_fraction << ")" << std::endl;
          if (norm_fraction < fractionCutoff_) {
	    auto clusterIndex_it = std::find(clusterIndex.begin(), clusterIndex.end(), clusterIndex[j]);
	    auto fraction_it = std::find(fractions.begin(), fractions.end(), fractions[j]);
	    clusterIndex.erase(clusterIndex_it);
	    fractions.erase(fraction_it);
	  }
	  /*if (norm_fraction >= fractionCutoff_){
            cellsIdInCluster[clusterIndex[j]].push_back(
              std::make_pair(i, norm_fraction));
          }*/
	}
	auto tot_norm_cleaned_fractions = std::accumulate(fractions.begin(), fractions.end(), 0.);
	for (unsigned int j = 0; j < clusterIndex.size(); j++) {
	  float norm_fraction = fractions[j]/tot_norm_cleaned_fractions;
	  cellsIdInCluster[clusterIndex[j]].push_back(
	    std::make_pair(i, norm_fraction));
	}  
      }
    }

    std::vector<std::pair<DetId, float>> thisCluster;
    for (int clIndex = 0; clIndex < numberOfClustersPerLayer_[layerId]; clIndex++){
      auto& cl = cellsIdInCluster[clIndex];
      auto position = calculatePosition(cl, layerId);
      float energy = 0.f;
      int seedDetId = -1;

      for (auto [cellIdx, fraction] : cl) {
	energy += cellsOnLayer.weight[cellIdx]*fraction;
	thisCluster.emplace_back(cellsOnLayer.detid[cellIdx], fraction);
	if (cellsOnLayer.isSeed[cellIdx]) {
	  seedDetId = cellsOnLayer.detid[cellIdx];
	}
      }
      auto globalClusterIndex = clIndex  + firstClusterIdx;

      if constexpr (std::is_same_v<T, EBLayerTiles>) {
	clusters_v_[globalClusterIndex] = 
	  reco::BasicCluster(energy, position, reco::CaloID::DET_ECAL_BARREL, thisCluster, algoId_);
      } else if constexpr (std::is_same_v<T, HBLayerTiles>) {
	  clusters_v_[globalClusterIndex] =
	    reco::BasicCluster(energy, position, reco::CaloID::DET_HCAL_BARREL, thisCluster, algoId_);
      } else {
	  clusters_v_[globalClusterIndex] = 
	    reco::BasicCluster(energy, position, reco::CaloID::DET_HO, thisCluster, algoId_);
      }
      clusters_v_[globalClusterIndex].setSeed(seedDetId);
      thisCluster.clear();
    }
    cellsIdInCluster.clear();
  }
  return clusters_v_;
}

template <typename T>
  math::XYZPoint SimBarrelCLUEAlgoT<T>::calculatePosition(const std::vector<std::pair<int, float>>& v, const unsigned int layerId) const {
  float total_weight = 0.f;
  float x = 0.f;
  float y = 0.f;
  float z = 0.f;

  auto& cellsOnLayer = cells_[layerId];
  for (const auto & [i, fraction] : v) {
    float rhEnergy = cellsOnLayer.weight[i] * fraction;
    total_weight += rhEnergy;
    float theta = 2 * std::atan(std::exp(-cellsOnLayer.eta[i]));
    const GlobalPoint gp(GlobalPoint::Polar(theta, cellsOnLayer.phi[i], cellsOnLayer.r[i]));
    x += gp.x() * rhEnergy;
    y += gp.y() * rhEnergy;
    z += gp.z() * rhEnergy;
  }

  if (total_weight != 0.f) {
    float inv_total_weight = 1.f/total_weight;
    return math::XYZPoint(x * inv_total_weight, y * inv_total_weight, z * inv_total_weight);
  } else {
    return math::XYZPoint(0.f, 0.f, 0.f);
  }
}

template <typename T>
void SimBarrelCLUEAlgoT<T>::calculateLocalDensity(const T& lt, const unsigned int layerId, float delta_c, float delta_r) {
  auto& cellsOnLayer = cells_[layerId];
  unsigned int numberOfCells = cellsOnLayer.detid.size();

  for (unsigned int i = 0; i < numberOfCells; ++i) {
    float delta = delta_c;
    std::array<int, 4> search_box = lt.searchBox(cellsOnLayer.eta[i] - delta,
						       cellsOnLayer.eta[i] + delta, 
						       cellsOnLayer.phi[i] - delta,
						       cellsOnLayer.phi[i] + delta);
    for (int etaBin = search_box[0]; etaBin < search_box[1]; ++etaBin) {
      for (int phiBin = search_box[2]; phiBin < search_box[3]; ++phiBin) {
	int phi = (phiBin % T::type::nRows);
	int binId = lt.getGlobalBinByBin(etaBin, phi);
	size_t binSize = lt[binId].size();

	for (unsigned int j = 0; j < binSize; ++j) {
	  unsigned int otherId = lt[binId][j];
	  if (distance(i, otherId, layerId) < delta) {
	    cellsOnLayer.rho[i] += (i == otherId ? 1.f : 0.5f) * cellsOnLayer.weight[otherId];
	  }
	}
      }
    }
  }
}

template <typename T>
void SimBarrelCLUEAlgoT<T>::calculateDistanceToHigher(const T& lt, const unsigned int layerId, float delta_c, float delta_r) {
  auto& cellsOnLayer = cells_[layerId];
  unsigned int numberOfCells = cellsOnLayer.detid.size();

  for (unsigned int i = 0; i < numberOfCells; ++i) {
    float maxDelta = std::numeric_limits<float>::max();
    float i_delta = maxDelta;
    float i_nearestHigher = -1;

    //float delta = delta_c;
    //auto range = outlierDeltaFactor_;
    auto range = delta_c;
    std::array<int, 4> search_box = lt.searchBox(cellsOnLayer.eta[i] - range,
						       cellsOnLayer.eta[i] + range,
						       cellsOnLayer.phi[i] - range,
						       cellsOnLayer.phi[i] + range);

    for (int etaBin = search_box[0]; etaBin < search_box[1]; ++etaBin) {
      for (int phiBin = search_box[2]; phiBin < search_box[3]; ++phiBin) {
	int phi = (phiBin % T::type::nRows);
	size_t binId = lt.getGlobalBinByBin(etaBin, phi);
	size_t binSize = lt[binId].size();
	
	for (unsigned int j = 0; j < binSize; ++j) {
	  int otherId = lt[binId][j];
	  float dist = distance(i, otherId, layerId);
	  bool foundHigher = (cellsOnLayer.rho[otherId] > cellsOnLayer.rho[i]) ||
			     (cellsOnLayer.rho[otherId] == cellsOnLayer.rho[i] &&
			     cellsOnLayer.detid[otherId] > cellsOnLayer.detid[i]);
	  if (foundHigher && dist <= i_delta) {
	    i_delta = dist;
	    i_nearestHigher = otherId;
	  }
	}
      }
    }
    bool foundNearestHigherInSearchBox = (i_delta != maxDelta);
    if (foundNearestHigherInSearchBox) {
      cellsOnLayer.delta[i] = i_delta;
      cellsOnLayer.nearestHigher[i] = i_nearestHigher;
    } else {
      cellsOnLayer.delta[i] = maxDelta;
      cellsOnLayer.nearestHigher[i] = -1;
    }
  }
}

template <typename T>
int SimBarrelCLUEAlgoT<T>::findAndAssignClusters(const unsigned int layerId, float delta_c, float delta_r) {
  unsigned int nClustersOnLayer = 0;
  auto& cellsOnLayer = cells_[layerId];
  unsigned int numberOfCells = cellsOnLayer.detid.size();
  std::vector<int> localStack;
  for (unsigned int i = 0; i < numberOfCells; ++i) {
    //float rho_c = rhoc_; //for testing purposes
    float rho_c = kappa_ * cellsOnLayer.sigmaNoise[i];
    float delta = delta_c;
    // cellsOnLayer.clusterIndex[i] = -1;
    bool isSeed = (cellsOnLayer.delta[i] > delta) && (cellsOnLayer.rho[i] >= rho_c);
    bool isOutlier = (cellsOnLayer.delta[i] > outlierDeltaFactor_ ) && (cellsOnLayer.rho[i] < rho_c);
    if (isSeed) {
      cellsOnLayer.clusterIndex[i].push_back(nClustersOnLayer);
      cellsOnLayer.isSeed[i] = true;
      cellsOnLayer.seedToCellIndex.push_back(i);
      nClustersOnLayer++;
      localStack.push_back(i);
    } else if (!isOutlier) {
      cellsOnLayer.followers[cellsOnLayer.nearestHigher[i]].push_back(i);
    }
  }

  while (!localStack.empty()) {
    int endStack = localStack.back();
    auto& thisSeed = cellsOnLayer.followers[endStack];
    localStack.pop_back();

    for (int j : thisSeed) {
      
      cellsOnLayer.clusterIndex[j].push_back(cellsOnLayer.clusterIndex[endStack][0]);
      localStack.push_back(j);
    }
  }
  return nClustersOnLayer;
}

template <typename T>
void SimBarrelCLUEAlgoT<T>::passSharedClusterIndex(const T& lt, const unsigned int layerId, float delta_c) {
  auto& cellsOnLayer = cells_[layerId];
  unsigned int numberOfCells = cellsOnLayer.detid.size();
  float delta = delta_c;
  for (unsigned int i = 0; i < numberOfCells; ++i) {
    // Do not run on outliers, but run also on seeds and followers
    if ((cellsOnLayer.clusterIndex[i].size() == 0)) continue;

    std::array<int, 4> search_box = lt.searchBox(cellsOnLayer.eta[i] - delta,
						       cellsOnLayer.eta[i] + delta,
						       cellsOnLayer.phi[i] - delta, 
						       cellsOnLayer.phi[i] + delta);
    
    for (int etaBin = search_box[0]; etaBin < search_box[1]; ++etaBin) {
      for (int phiBin = search_box[2]; phiBin < search_box[3]; ++phiBin) {
	int phi = (phiBin % T::type::nRows);
	size_t binId = lt.getGlobalBinByBin(etaBin, phi);
	size_t binSize = lt[binId].size();

	for (unsigned int j = 0; j < binSize; ++j) {
	  unsigned int otherId = lt[binId][j];
	  if (cellsOnLayer.clusterIndex[otherId].size() == 0) continue;
	  int otherClusterIndex = cellsOnLayer.clusterIndex[otherId][0];
	  if (std::find(std::begin(cellsOnLayer.clusterIndex[i]),
                        std::end(cellsOnLayer.clusterIndex[i]), otherClusterIndex) == std::end(cellsOnLayer.clusterIndex[i]))
            cellsOnLayer.clusterIndex[i].push_back(otherClusterIndex);
	}	  
      }
    }
  }
}

template <typename T> 
void SimBarrelCLUEAlgoT<T>::setDensity(const unsigned int layerId) {
  auto& cellsOnLayer = cells_[layerId];
  unsigned int numberOfCells = cellsOnLayer.detid.size();
  for (unsigned int i = 0; i < numberOfCells; ++i) {
    density_[cellsOnLayer.detid[i]] = cellsOnLayer.rho[i];
  }
}

template <typename T>
Density SimBarrelCLUEAlgoT<T>::getDensity() {
  return density_;
}

template class SimBarrelCLUEAlgoT<EBLayerTiles>;
template class SimBarrelCLUEAlgoT<HBLayerTiles>;
template class SimBarrelCLUEAlgoT<HOLayerTiles>;

