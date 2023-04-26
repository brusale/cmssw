#ifndef RecoLocalCalo_HGCalRecProducers_BarrelCLUEAlgo_h
#define RecoLocalCalo_HGCalRecProducers_BarrelCLUEAlgo_h

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalClusteringAlgoBase.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalLayerTiles.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "CondFormats/DataRecord/interface/EcalPFRecHitThresholdsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalPFRecHitThresholds.h"
#include "CondFormats/DataRecord/interface/HcalPFCutsRcd.h"
#include "CondFormats/HcalObjects/interface/HcalPFCuts.h"

// C/C++ headers
#include <set>
#include <string>
#include <vector>

using Density = hgcal_clustering::Density;

template <typename TILE>
class BarrelCLUEAlgoT : public HGCalClusteringAlgoBase {
public:
  BarrelCLUEAlgoT(const edm::ParameterSet& ps, edm::ConsumesCollector iC)
      : HGCalClusteringAlgoBase(
            (HGCalClusteringAlgoBase::VerbosityLevel)ps.getUntrackedParameter<unsigned int>("verbosity", 3),
            reco::CaloCluster::undefined,
            iC),
        tok_ebThresholds_(iC.esConsumes<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd>()),
	tok_hcalThresholds_(iC.esConsumes<HcalPFCuts, HcalPFCutsRcd>()),
	vecDeltas_(ps.getParameter<std::vector<double>>("deltac")),
	rhoc_(ps.getParameter<double>("rhoc")),
	maxLayerIndex_(ps.getParameter<int>("maxLayerIndex")) {}
  ~BarrelCLUEAlgoT() override {}

  void getEventSetupPerAlgorithm(const edm::EventSetup& es) override;

  void populate(const HGCRecHitCollection& hits) override {};
  void populate(const reco::PFRecHitCollection& hits) override;

  // this is the method that will start the clusterisation (it is possible to invoke this method
  // more than once - but make sure it is with different hit collections (or else use reset)

  void makeClusters() override;

  // this is the method to get the cluster collection out
  std::vector<reco::BasicCluster> getClusters(bool) override;

  void reset() override {
    clusters_v_.clear();
    clusters_v_.shrink_to_fit();
    for (auto& cl : numberOfClustersPerLayer_) {
      cl = 0;
    }

    for (auto& cells : cells_) {
      cells.clear();
      cells.shrink_to_fit();
    }
    density_.clear();
  }

  Density getDensity() override;

  void computeThreshold();

  static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {
    iDesc.add<int>("maxLayerIndex");
    iDesc.add<double>("rhoc");
    iDesc.add<std::vector<double>>("deltac",
                                   {
                                       0.0175,
                                       5*0.087,
                                       5*0.087
                                   });
  }

  /// point in the space
  typedef math::XYZPoint Point;

private:
  // To get ECAL sigmaNoise
  edm::ESGetToken<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd> tok_ebThresholds_;
  const EcalPFRecHitThresholds* ebThresholds_;
  //To get HCAL sigmaNoise
  edm::ESGetToken<HcalPFCuts, HcalPFCutsRcd> tok_hcalThresholds_;
  const HcalPFCuts* hcalThresholds_;
  // The two parameters used to identify clusters
  std::vector<double> vecDeltas_;
  double kappa_;
  double rhoc_;
  int maxLayerIndex_;

  Density density_;
  // For keeping the density per hit


  float outlierDeltaFactor_ = 2.f;

  struct BarrelCellsOnLayer {
    std::vector<DetId> detid;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> r;

    std::vector<float> weight;
    std::vector<float> rho;

    std::vector<float> delta;
    std::vector<int> nearestHigher;
    std::vector<int> clusterIndex;
    std::vector<float> sigmaNoise;
    std::vector<std::vector<int>> followers;
    std::vector<bool> isSeed;

    void clear() {
      detid.clear();
      eta.clear();
      phi.clear();
      r.clear();
      weight.clear();
      rho.clear();
      delta.clear();
      nearestHigher.clear();
      clusterIndex.clear();
      sigmaNoise.clear();
      followers.clear();
      isSeed.clear();
    }

    void shrink_to_fit() {
      detid.shrink_to_fit();
      r.shrink_to_fit();
      eta.shrink_to_fit();
      phi.shrink_to_fit();
      weight.shrink_to_fit();
      rho.shrink_to_fit();
      delta.shrink_to_fit();
      nearestHigher.shrink_to_fit();
      clusterIndex.shrink_to_fit();
      sigmaNoise.shrink_to_fit();
      followers.shrink_to_fit();
      isSeed.shrink_to_fit();
    }
  };

  std::vector<BarrelCellsOnLayer> cells_;

  std::vector<int> numberOfClustersPerLayer_;

  inline float distance2(int cell1, int cell2, int layerId) const {  // distance squared
    const float dphi = reco::deltaPhi(cells_[layerId].phi[cell1], cells_[layerId].phi[cell2]);
    const float deta = cells_[layerId].eta[cell1] - cells_[layerId].eta[cell2];
    return (deta * deta + dphi * dphi);
  }

  inline float distance(int cell1, int cell2, int layerId) const {  // 2-d distance on the layer (x-y)
    return std::sqrt(distance2(cell1, cell2, layerId));
  }

  void prepareDataStructures(const unsigned int layerId);
  void calculateLocalDensity(const TILE& lt,
                             const unsigned int layerId,
                             float delta_c,
                             float delta_r);  // return max density
  void calculateDistanceToHigher(const TILE& lt, const unsigned int layerId, float delta_c, float delta_r);
  int findAndAssignClusters(const unsigned int layerId, float delta_c, float delta_r);
  math::XYZPoint calculatePosition(const std::vector<int>& v, const unsigned int layerId) const;
  void setDensity(const unsigned int layerId);
};

// explicit template instantiation
extern template class BarrelCLUEAlgoT<EBLayerTiles>;
extern template class BarrelCLUEAlgoT<HBLayerTiles>;
extern template class BarrelCLUEAlgoT<HOLayerTiles>;

using EBCLUEAlgo = BarrelCLUEAlgoT<EBLayerTiles>;
using HBCLUEAlgo = BarrelCLUEAlgoT<HBLayerTiles>;
using HOCLUEAlgo = BarrelCLUEAlgoT<HOLayerTiles>;

#endif
