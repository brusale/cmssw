#include "RecoLocalCalo/HGCalRecProducers/plugins/BarrelCLUEAlgo.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"


#include "oneapi/tbb/task_arena.h"
#include "oneapi/tbb.h"
#include <limits>

using namespace hgcal_clustering;

template <typename T>
void BarrelCLUEAlgoT<T>::getEventSetupPerAlgorithm(const edm::EventSetup& es) {
  cells_.clear();
  numberOfClustersPerLayer_.clear();
  cells_.resize(maxlayer_ + 1);
  numberOfClustersPerLayer_.resize(maxlayer_ + 1, 0);
}

template <typename T>
void BarrelCLUEAlgoT<T>::populate(const reco::PFRecHitCollection& hits) {
  for (unsigned int i = 0; i < hits.size(); ++i) {
    const reco::PFRecHit hit = hits[i];
    DetId detid(hit.detId());
    const GlobalPoint position = rhtools_.getPosition(detid);
    int layer = 0;

    if (detid.det() == DetId::Hcal) {
      HcalDetId hid(detid);
      layer = hid.depth();
      if (detid.subdetId() == HcalSubdetector::HcalOuter) 
	layer += 1;
    }

    cells_[layer].detid.emplace_back(detid);
    cells_[layer].eta.emplace_back(position.eta());
    cells_[layer].phi.emplace_back(position.phi());
    cells_[layer].weight.emplace_back(hit.energy());
    cells_[layer].r.emplace_back(position.mag());
    //cells_[layer].sigmaNoise.emplace_back(sigmaNoise);
  }
}

template <typename T>
void BarrelCLUEAlgoT<T>::prepareDataStructures(unsigned int l) {
  auto cellsSize = cells_[l].detid.size();
  cells_[l].rho.resize(cellsSize, 0.f);
  cells_[l].delta.resize(cellsSize, 9999999);
  cells_[l].nearestHigher(cellsSize, -1);
  cells_[l].clusterIndex.resize(cellsSize, -1);
  cells_[l].followers.resize(cellsSize);
  cells_[l].eta.resize(cellsSize, 0.f);
  cells_[l].phi.resize(cellsSize, 0.f);
  cells_[l].r.resize(cellsSize, 0.f);
}

template <typename T>
void BarrelCLUEAlgoT<T>::makeClusters() {
  tbb::this_task_arena::isolate([&] {
    tbb::parallel_for(size_t(0), size_t(maxlayer_ + 1), [&](size_t i) {
      prepareDataStructures(i);
      T lt;
      lt.clear();
      lt.fill(cells_[i].eta, cells_[i].phi);
      float delta_c;
      
      if (i == 0)
	delta_c = vecDeltas_[0];
      else
	delta_c = vecDeltas_[1];
      float delta_r = vecDeltas_[2];
      
      calculateLocalDensity(lt, i, delta_c, delta_r);
      calculateDistanceToHigher(lt, i, delta_c, delta_r);
      numberOfClustersPerLayer_[i] = findAndAssignClusters(i, delta_c, delta_r);
      });
    });
  for (unsigned int i = 0; i < maxlayer_ + 1; ++i) {
    setDensity(i);
  }
}

template <typename T>
std::vector<reco::BasicCluster> BarrelCLUEAlgoT<T>::getClusters(bool) {
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
  std::vector<std::vector<int>> cellsIdInCluster;
  cellsIdInCluster.reserve(maxClustersOnLayer);

  for (unsigned int layerId = 0; layerId < maxlayer_ + 1; ++layerId) {
    cellsIdInCluster.resize(numberOfClustersPerLayer_[layerId]);
    auto& cellsOnLayer = cells_[layerId];
    unsigned int numberOfCells = cellsOnLayer.detid.size();
    auto firstClusterIdx = offsets[layerId];

    for (unsigned int i = 0; i < numberOfCells; ++i) {
      auto clusterIndex = cellsOnLayer.clusterIndex[i];
      if (clusterIndex != -1) 
	cellsIdInCluster[clusterIndex].push_back(i);
    }

    std::vector<std::pair<DetId, float>> thisCluster;

    for (auto& cl : cellsIdInCluster) {
      auto position = calculatePosition(cl, layerId);
      float energy = 0.f;
      int seedDetId = -1;

      for (auto cellIdx : cl) {
	energy += cellsOnLayer.weight[cellIdx];
	thisCluster.emplace_back(cellsOnLayer.detid[cellIdx], 1.f);
	if (cellsOnLayer.isSeed[cellIdx]) {
	  seedDetId = cellsOnLayer.detid[cellIdx];
	}
      }
      auto globalClusterIndex = cellsOnLayer.clusterIndex[cl[0]] + firstClusterIdx;

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
math::XYZPoint BarrelCLUEAlgoT<T>::calculatePosition(const std::vector<int>& v, const unsigned int layerId) const {
  float total_weight = 0.f;
  float x = 0.f;
  float y = 0.f;
  float z = 0.f;

  auto& cellsOnLayer = cells_[layerId];
  for (auto i : v) {
    float rhEnergy = cellsOnLayer.weight[i];
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
