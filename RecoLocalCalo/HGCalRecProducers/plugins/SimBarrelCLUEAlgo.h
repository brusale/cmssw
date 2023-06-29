#ifndef RecoLocalCalo_HGCalRecProducers_SimBarrelCLUEAlgo_h
#define RecoLocalCalo_HGCalRecProducers_SimBarrelCLUEAlgo_h

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

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
// C/C++ headers
#include <set>
#include <string>
#include <vector>

using Density = hgcal_clustering::Density;

template <typename TILE>
class SimBarrelCLUEAlgoT : public HGCalClusteringAlgoBase {
public:
  SimBarrelCLUEAlgoT(const edm::ParameterSet& ps, edm::ConsumesCollector iC)
      : HGCalClusteringAlgoBase(
            (HGCalClusteringAlgoBase::VerbosityLevel)ps.getUntrackedParameter<unsigned int>("verbosity", 3),
            reco::CaloCluster::undefined),
        tok_ebThresholds_(iC.esConsumes<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd>()),
	tok_hcalThresholds_(iC.esConsumes<HcalPFCuts, HcalPFCutsRcd>()),
	vecDeltas_(ps.getParameter<std::vector<double>>("deltac")),
	kappa_(ps.getParameter<double>("kappa")),
	rhoc_(ps.getParameter<double>("rhoc")),
        fractionCutoff_(ps.getParameter<double>("fractionCutoff")),
	maxLayerIndex_(ps.getParameter<int>("maxLayerIndex")),
	outlierDeltaFactor_(ps.getParameter<double>("outlierDeltaFactor")) {}
  ~SimBarrelCLUEAlgoT() override {}

  void getEventSetupPerAlgorithm(const edm::EventSetup& es) override;

  void populate(const HGCRecHitCollection& hits) override {};
  void populate(const reco::PFRecHitCollection& hits) override {};
  void populate(const edm::PCaloHitContainer& hits) override;
  //void populate(const std::vector<std::pair<uint32_t, float>> &hae);

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
    iDesc.add<double>("outlierDeltaFactor", 2.);
    iDesc.add<double>("kappa", 1.34);
    iDesc.add<int>("maxLayerIndex");
    iDesc.add<double>("rhoc");
    iDesc.add<std::vector<double>>("deltac",
                                   {
                                       0.0175,
                                       5*0.087,
                                       5*0.087
                                   });
    iDesc.add<double>("fractionCutoff", 0.);
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
  double fractionCutoff_;
  int maxLayerIndex_;

  Density density_;
  // For keeping the density per hit


  float outlierDeltaFactor_;
  struct BarrelCellsOnLayer {
    std::vector<DetId> detid;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> r;

    std::vector<float> weight;
    std::vector<float> rho;

    std::vector<float> delta;
    std::vector<int> nearestHigher;
    //std::vector<int> clusterIndex;
    std::vector<std::vector<int>> clusterIndex;
    std::vector<float> sigmaNoise;
    std::vector<std::vector<int>> followers;
    std::vector<bool> isSeed;
    std::vector<int> seedToCellIndex;

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
      seedToCellIndex.clear();
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
      seedToCellIndex.shrink_to_fit();
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
  void passSharedClusterIndex(const TILE& lt, const unsigned int layerId, float delta_c);
  math::XYZPoint calculatePosition(const std::vector<std::pair<int,float>>& v, const unsigned int layerId) const;
  void setDensity(const unsigned int layerId);
};

// explicit template instantiation
extern template class SimBarrelCLUEAlgoT<EBLayerTiles>;
extern template class SimBarrelCLUEAlgoT<HBLayerTiles>;
extern template class SimBarrelCLUEAlgoT<HOLayerTiles>;

using SimEBCLUEAlgo = SimBarrelCLUEAlgoT<EBLayerTiles>;
using SimHBCLUEAlgo = SimBarrelCLUEAlgoT<HBLayerTiles>;
using SimHOCLUEAlgo = SimBarrelCLUEAlgoT<HOLayerTiles>;

#endif
