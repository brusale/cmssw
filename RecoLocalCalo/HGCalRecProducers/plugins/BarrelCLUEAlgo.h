#ifndef RecoLocalCalo_HGCalRecProducers_BarrelCLUEAlgo_h
#define RecoLocalCalo_HGCalRecProducers_BarrelCLUEAlgo_h

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalClusteringAlgoBase.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/EBTilesConstants.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HBTilesConstants.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HOTilesConstants.h"

using Density = hgcal_clustering::Density;

template <typename TILE>
class BarrelCLUEAlgoT : HGCalClusteringAlgoBase {
  public:
    BarrelCLUEAlgoT(const edm::ParameterSet& ps, edm::ConsumesCollector iC) 
      : HGCalClusteringAlgoBase(
	          (HGCalClusteringAlgoBase::VerbosityLevel)ps.getUntrackedParameter<unsigned int>("verbosity", 3);
	          reco::CaloCluster::undeined,
	          iC) {}
    ~BarrelCLUEAlgoT() override {}

    void getEventSetupPerAlgorithm(const edm::EventSetup& es) override;

    void populate(reco::PFRecHitCollection& hits) override;

    void makeClusters() override;

    std::vector<reco::BasicCluster> getClusters(bool) override;

    void reset() override {
      clusters_v_.clear();
      clusters_v.shrink_to_fit();

      for (auto& cl : numberOfClustersPerLayer_) {
	      cl = 0;
      }

      for (auto& cells : cells_) {
	      cells.clear();
	      cells.shrink_to_fit();
      }
      density.clear();
    }


  private:
    std::vector<double> vecDeltas_;
    double kappa_;

    Density density_;

    float outlierDeltaFactor_ = 2.f;

    struct BarrelCellsOnLayer {
      std::vector<DetId> detid;
      //std::vector<float> eta;
      std::vector<float> phi;
      std::vector<float> x;
      std::vector<float> y;
      std::vector<float> z;
      std::vector<float> r;

      std::vector<float> weight;
      std::vector<float> rho;
      std::vector<float> delta;
      std::vector<int> nearestHigher;
      std::vector<int> clusterIndex;
      std::vector<std::vector<int>> followers;
      std::vector<bool> isSeed;
      //std::vector<float> sigmaNoise;

      void clear() {
        detid.clear();
        //eta.clear();
        phi.clear();
        x.clear();
        y.clear();
        z.clear();
        r.clear();
        
        weight.clear();
        rho.clear();
        delta.clear();
        nearestHigher.clear();
        clusterIndex.clear();
        followers.clear();
        isSeed.clear();
        //sigmaNoise.clear();
       }

      void shrink_to_fit() {
        detid.shrink_to_fit();
        //eta.shrink_to_fit();
        phi.shrink_to_fit();
        z.shrink_to_fit();
        r.shrink_to_fit();
        
        weight.shrink_to_fit();
        rho.shrink_to_fit();
        delta.shrink_to_fit();
        nearestHigher.shrink_to_fit();
        clusterIndex.shrink_to_fit();
        followers.shrink_to_fit();
        isSeed.shrink_to_fit();
        //sigmaNoise.shrink_to_fit();
       }
    };

    struct HcalCellsOnLayer : BarrelCellsOnLayer {
      std::vector<int> depth;

      void clear() {
	      BarrelCellsOnLayer::clear();
	      depth.clear();
      }

      void shrink_to_fit() {
	      BarrelCellsOnLayer::shrink_to_fit();
	      depth.shrink_to_fit();
      }
    };

    if constexpr (std::is_same_v<TILE, EBLayerTiles>) {
      std::vector<BarrelCellsOnLayer> cells_;
    } else {
      std::vector<HcalCellsOnLayer> cells_;
    }

    std::vector<int> numberOfClustersPerLayer_;

    inline float distance2(int cell1, int cell2, int layerId) {
      //const float dphi = reco::deltaPhi(cells_[layerId].phi[cell1] - cells_[layerId].phi[cell2]);
      //const float deta = cells_[layerId].eta[cell1] - cells_[layerId].eta[cell2];
      const float drphi = r * reco::deltaPhi(cells_[layerId].phi[cell1] - cells_[layerId].phi[cell2]);
      const float dz = cells_[layerId].z[cell1] - cells_[layerId].z[cell2];
      //return (deta * deta + dphi * dphi)
      return (drphi * drphi + dz * dz);
    }

    inline float distance(int cell1, int cell2, int layerId) {
      return std::sqrt(distance2(cell1, cell2, layerId));
    }

    void prepareDataStructures(const unsigned int layerId);
    void calculateLocalDensity(const TILE& lt,
                               const unsigned int layerId,
                               float delta_c, 
                               float delta_r);
    void calculateDistanceToHigher(const TILE& lt, const unsigned int layerId, float delta_c, float delta_r);
    int findAndAssignClusters(const unsigned int layerId, float delta_c, float delta_r);
    math::XYZPoint calculatePosition(const std::vector<int>& v, const unsigned int layerId) const;
    void setDensity(const unsigned int layerId);

};

#endif
