#include "RecoLocalCalo/HGCalRecProducers/plugins/BarrelCLUEAlgo.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "oneapi/tbb/task_arena.h"
#include "oneapi/tbb.h"
#include <limits>

using namespace hgcal_clustering;

template <typename T>
void BarrelCLUEAlgoT<T>::getEventSetupPerAlgorithm(const edm::EventSetup& es) {
  cells_.clear();
  numberOfClustersPerLayer_.clear();
  cells_.resize(maxlayer_ + 1); //ECAL maxlayer = 0; HCAL maxlayer=6?
  numberOfClustersPerLayer_.resize(maxlayer_ + 1, 0);
}

template <typename T>
void BarrelCLUEAlgoT<T>::populate(const reco::PFRecHitCollection& hits) {
  for (unsigned int i = 0; i < hits.size(); ++i) {
    const reco::PFRecHit& hit = hits[i];
    DetId detid(hit.detId());
    const GlobalPoint position(rhtools_.getPosition(detid));
    int layer = 0;

    cells_[layer].detid.emplace_back(detid);
    //cells_[layer].eta.emplace_back(position.eta());
    cells_[layer].phi.emplace_back(position.phi());
    cells_[layer].r.emplace_back(position.mag());
    cells_[layer].x.emplace_back(position.x());
    cells_[layer].y.emplace_back(position.y());
    cells_[layer].z.emplace_back(position.z());
    
    cells_[layer].weight.emplace_back(hit.energy());
    if (detid.det() == DetId::Hcal) {
      HcalDetId hid(hit.detid());
      layer = hid.depth();
      if (detid.subdetId() == HcalSubdetector::HcalOuter) {
	      layer += 1;
      }
      cells_[layer].depth.emplace_back(layer);
    }
  }
}

template <typename T>
void BarrelCLUEAlgoT<T>::prepareDataStructures(unsigned int l) {
  auto cellsSize = cells_[l].detid.size();
  cells_[l].rho.resize(cellsSize, 0.f);
  cells_[l].delta.resize(cellsSize, 9999999);
  cells_[l].nearestHigher.resize(cellsSize, -1);
  cells_[l].clusterIndex.resize(cellsSize, -1);
  cells_[l].followers.resize(cellsSize);
  cells_[l].isSeed.resize(cellsSize, false);
  //cells_[l].eta.resize(cellsSize, true);
  cells_[l].phi.resize(cellsSize, 0.f);
  cells_[l].x.resize(cellsSize, 0.f);
  cells_[l].y.resize(cellsSize, 0.f);
  cells_[l].z.resize(cellsSize, 0.f); 
  cells_[l].r.resize(cellsSize, 0.f);

}

template <typename T>
void BarrelCLUEAlgoT<T>::makeClusters() {
  tbb::this_task_arena::isolate([&] {
    tbb::parallel_for(size_t(0), size_t(maxlayer_ + 1), [&](size_t i) {
      prepareDataStructures(i);
      T lt;
      lt.clear();
      //lt.fill(cells_[i].eta, cells_[i].phi)
      //lt.fill(cells_[i].z, cells_[i].phi, cells_[i].x, cells_[i].y);
      lt.fill(cells_[i].z, cells_[i].phi);
      float delta_c;
      if (i == 0) {  //need to check this
        delta_c = vecDeltas_[0];
      } else {
        delta_c = vecDeltas_[1];
      }
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

