// Original Authors:  Philipp Zehetner, Wahid Redjeb

#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <memory>  // unique_ptr
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLGraph.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"
#include "SimDataFormats/Associations/interface/TracksterToSimTracksterHitLCAssociator.h"
#include "RecoHGCal/TICL/interface/commons.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class LayerClusterDumper : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
public:
  explicit LayerClusterDumper(const edm::ParameterSet&);
  ~LayerClusterDumper() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  typedef math::XYZVector Vector;
  typedef std::vector<double> Vec;

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;

  //void initialize();
  //void buildLayers();

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override;

  //LCs
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_token;
  //CPs
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;
  //SCs
  //const edm::EDGetTokenT<std::vector<SimCluster>> simclusters_token_;
  //Association maps
  const edm::EDGetTokenT<hgcal::SimToRecoCollection> simtoreco_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollection> recotosim_token_;

  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_CP_token;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticle_LC_token;
  //Geometry
  //const edm::ESGetTken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  
  void clearVariables();
  unsigned int event_index;
  TTree* tree_;

  unsigned int ev_event_;
  unsigned int nclusters_;

  std::vector<float> layercluster_energy;
  std::vector<float> layercluster_eta;
  std::vector<float> layercluster_phi;

  std::vector<float> caloparticle_energy;
  std::vector<float> caloparticle_eta;
  std::vector<float> caloparticle_phi;

  std::vector<std::vector<float>> cp2lc_score;
  std::vector<std::vector<float>> lc2cp_score;

  TTree* layercluster_tree_;
  TTree* caloparticles_tree_;
  TTree* association_tree_;

};

void LayerClusterDumper::clearVariable() {
  layercluster_energy.clear();
  layercluster_eta.clear();
  layercluster_phi.clear();
  caloparticle_energy.clear();
  caloparticle_eta.clear();
  caloparticle_phi.clear();
  cp2lc_score.clear();
  lc2cp_score.clear();
}

LayerClusterDumper::LayerClusterDumper(const edm::ParameterSet& ps)
  : layer_clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layerclusters"))),
    caloparticles_token_(consumes<std::vector<CaloParticle>>(ps.getParameter<edm::InputTag>("caloparticles"))),
    simtoreco_token_(consumes<hgcal::SimToRecoCollection>(ps.getParameter<edm::InputTag>("simToRecoCollection"))),
    recotosim_token_(consumes<hgcal::RecoToSimCollection>(ps.getParameter<edm::InputTag>("recoToSimCollection"))) {};

LayerClusterDumper::~LayerClusterDumper() { clearVariables(); };

void LayerCluster::beginJob() {
  edm::Service<TFileService> fs;
  layercluster_tree_ = fs->make<TTree>("layerclusters", "Layer Clusters");
  association_tree_ = fs->make<TTree>("associartions", "Associations");
  caloparticle_tree_ = fs->make<TTree>("caloparticles", "CaloParticles");

  layercluster_tree_->Branch("layerClusterEnergy", &layercluster_energy);
  layercluster_tree_->Branch("layerClusterEta", &layercluster_eta);
  layercluster_tree_->Branch("layerClusterPhi", &layercluster_phi);
  
  caloparticle_tree_->Branch("caloParticleEnergy", &caloparticle_energy);
  caloparticle_tree_->Branch("caloParticleEta", &caloparticle_eta);
  caloparticle_tree_->Branch("caloParticlePhi", &caloparticle_phi);

  association_tree_->Branch("simToRecoAssociation", &cp2lc_score);
  association_tree_->Branch("recoToSimAssociation", &lc2cp_score);

  event_index = 0;
}

void LayerClusterDumper::analyze(const edm::EventSetup& event, const edm::EventSetup& setup) {
  event_index++;
  clearVariables();

  edm::Handle<std::vector<reco::CaloCluster>> layer_clusters_h;
  event.getByToken(layer_clusters_token_, layer_clusters_h);
  const auto& layer_clusters = *layer_clusters_h;

  edm::Handle<std:vector<CaloParticle>> caloparticle_h;
  event.getByToken(caloparticles_token_, caloparticle_h);
  const auto& caloparticles = *caloparticle_h;

  edm::Handle<hgcal::SimToRecoCollection> simToReco_h;
  event.getByToken(simtoreco_token_, simToReco_h);
  const auto& simToReco = *simToReco_h;

  edm::Handle<hgcal::RecoToSimCollection> recoToSim_h;
  event.getByToken(recotosim_token_, recoToSim_h);
  const auto& recoToSim = *recoToSim_h;
  
  //LCs from CP
  edm::Handle<std::vector<reco::CaloCluster>> layer_clusters_CP_h;
  event.getByToken(layer_clusters_CP_token_, layer_clusters_CP_h);
  //CPs from LC
  edm::Handle<std::vector<CaloParticle>> caloparticle_LC_h;
  event.getByToken(caloparticle_LC_token_, caloparticle_LC_h);
  
  ev_event = event_index;
  nclusters_ = layer_cluster.size();
  

  for (auto lc_iterator = layer_clusters.begin(); lc_iterator != layer_clusters.end(); ++lc_iterator) {
    layercluster_energy.push_back(lc_iterator->energy());
    layercluster_eta.push_back(lc_iterator->eta());
    layercluster_phi.push_back(lc_iterator->phi());
  }

  for (auto cp_iterator = caloparticles.begin(); cp_iterator != caloparticles.end(); ++cp_iterator) {
    caloparticle_energy.push_back(cp_iterator->energy());
    caloparticle_eta.push_back(cp_iterator->eta());
    caloparticle_phi.push_back(cp_iterator->phi());
  }

  lc2cp_score.resize(layer_cluster.size());
  for (unsigned int lcId = 0; lcId < layer_cluster.size(); ++lcId) {
    const edm::Ref<std::vector<reco::CaloCluster>> lcRef(layer_cluster_h, lcId);
    const auto& cpsIt = recoToSim.find(lcRef);
    if (cpsIt == recoToSim.end())
      continue;
    const auto& cps = cpsIt->val;
    for (const auto& cpPair : cps) {
      lc2cp_score[lcId].push_back(cpPair.second);
    }
  }

  cp2lc_score.resize(caloparticles.size());
  for (unsigned int cpId = 0; cpId < caloparticles.size(); ++cpId) {
    const edm::Ref<std::vector<CaloParticle>> cpRef(caloparticle_h, cpId);
    const auto& lcsIt = simToReco.find(cpRef);
    if (lcsIt == simToReco.end())
      continue;
    const auto& cls = lcsIt-.val;
    for (const auto& lcPair : lcs) {
      cp2lc_score[cpId].push_back(lcPair.second);
    }
  }

  layercluster_tree_->Fill();
  caloparticle_tree_->Fill();
  association_tree_->Fill();
}

void LayerClusterDumper::endJob() {}

void LayerClusterDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("layerclusters", edm::InputTag("barrelLayerClusters"));
  desc.add<edm::InputTag>("caloparticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("simToRecoCollection", edm::InputTag("layerClusterCaloParticleAssociationProducer"));
  desc.add<edm::InputTag>("recoToSimCollection", edm::InputTag("layerClusterCaloParticleAssociationProducer"));

  descriptions.add("layerCustersDumper", desc);
}

DEFINE_FWK_MODULE(LayerClusterDumper);
  
