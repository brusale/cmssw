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
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
//#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
//#include "TrackingTools/GeomPropagators/interface/Propagator.h"
//#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"
#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociatorBaseImpl.h"
#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociatorBaseImpl.h"
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

  const edm::EDGetTokenT<reco::PFRecHitCollection> pfrechit_token_;
  //LCs
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_token_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> sim_layer_clusters_token_;
  //CPs
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;
  //SCs
  const edm::EDGetTokenT<std::vector<SimCluster>> simclusters_token_;
  //Association maps
  const edm::EDGetTokenT<hgcal::SimToRecoCollection> simtoreco_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollection> recotosim_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollection> simtorecosim_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollection> recotosimsim_token_;

  const edm::EDGetTokenT<hgcal::SimToRecoCollectionWithSimClusters> simtoreco_simcl_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionWithSimClusters> recotosim_simcl_token_;
  //Geometry
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  
  void clearVariables();

  std::vector<std::pair<std::pair<float, float>, float>> addHitsAndDensities(std::vector<std::pair<DetId, float>> hits_and_energies, const CaloGeometry* geom, float totalEnergy);
  std::vector<std::pair<std::pair<float, float>, float>> addHitsAndDensities(std::vector<std::pair<uint32_t, float>> hits_and_fractions, const CaloGeometry* geom, float totalEnergy);

  float distance(GlobalPoint point1, GlobalPoint point2);
  unsigned int event_index;
  TTree* tree_;
  
  const CaloGeometry* geom;

  unsigned int ev_event_;
  unsigned int nclusters_;
  unsigned int ncaloparticles_;
  unsigned int nsimlayerclusters_;
  unsigned int npfrechits_;
  unsigned int nsimclusters_;

  std::vector<float> pfrechit_eta;
  std::vector<float> pfrechit_phi;
  std::vector<float> pfrechit_energy;
  std::vector<float> pfrechit_time;

  std::vector<float> layercluster_energy;
  std::vector<float> layercluster_eta;
  std::vector<float> layercluster_phi;
  std::vector<float> layercluster_seed_eta;
  std::vector<float> layercluster_seed_phi;
  std::vector<std::vector<float>> layercluster_hit_eta;
  std::vector<std::vector<float>> layercluster_hit_phi;
  std::vector<std::vector<float>> layercluster_hit_density;
  std::vector<std::vector<float>> layercluster_hit_energy;
  std::vector<int> layercluster_nhits;
  std::vector<std::vector<float>> seed_density;
  std::vector<float> layercluster_resolution;
  std::vector<int> layercluster_layer;
  std::vector<float> layercluster_pu_contribution;

  std::vector<float> simlayercluster_energy;
  std::vector<float> simlayercluster_eta;
  std::vector<float> simlayercluster_phi;
  std::vector<float> simlayercluster_seed_eta;
  std::vector<float> simlayercluster_seed_phi;
  std::vector<std::vector<float>> simlayercluster_hit_eta;
  std::vector<std::vector<float>> simlayercluster_hit_phi;
  std::vector<std::vector<float>> simlayercluster_hit_energy;
  std::vector<int> simlayercluster_nhits;
  std::vector<float> simlayercluster_resolution;

  std::vector<float> caloparticle_energy;
  std::vector<float> caloparticle_eta;
  std::vector<float> caloparticle_phi;
  std::vector<std::vector<float>> caloparticle_hit_eta;
  std::vector<std::vector<float>> caloparticle_hit_phi;
  std::vector<std::vector<float>> caloparticle_hit_density;
  std::vector<std::vector<float>> caloparticle_hit_energy;
  std::vector<std::vector<float>> shared_energy;
  std::vector<int> caloparticle_bx;

  std::vector<float> simcl_energy;
  std::vector<float> simcl_eta;
  std::vector<float> simcl_phi;
  std::vector<std::vector<float>> simcl_shared_energy;
  std::vector<int> simcl_layer;


  std::vector<int> layercluster_event;
  std::vector<int> layercluster_hit_event;

  std::vector<int> caloparticle_event;
  std::vector<int> caloparticle_hit_event;


  std::vector<std::vector<float>> cp2lc_score;
  std::vector<std::vector<uint32_t>> associatedcp_to_lc;
  std::vector<std::vector<float>> lc2cp_score;
  std::vector<std::vector<uint32_t>> associatedlc_to_cp;
  std::vector<std::vector<float>> cp2lc_score_sim;
  std::vector<std::vector<uint32_t>> associatedcp_to_simlc;
  std::vector<std::vector<float>> lc2cp_score_sim;
  std::vector<std::vector<uint32_t>> associatedsimlc_to_cp;
 
  std::vector<std::vector<float>> sc2lc_score;
  std::vector<std::vector<uint32_t>> associatedsc_to_lc;
  std::vector<std::vector<float>> lc2sc_score;
  std::vector<std::vector<uint32_t>> associatedlc_to_sc;

  TTree* layercluster_tree_;
  TTree* simlayercluster_tree_;
  TTree* caloparticle_tree_;
  TTree* simcluster_tree_;
  TTree* pfrechit_tree_;
};

void LayerClusterDumper::clearVariables() {
  layercluster_energy.clear();
  layercluster_eta.clear();
  layercluster_phi.clear();
  layercluster_seed_eta.clear();
  layercluster_seed_phi.clear();
  layercluster_layer.clear();
  layercluster_pu_contribution.clear();
  caloparticle_energy.clear();
  caloparticle_eta.clear();
  caloparticle_phi.clear();
  cp2lc_score.clear();
  lc2cp_score.clear();
  cp2lc_score_sim.clear();
  lc2cp_score_sim.clear();
  associatedcp_to_lc.clear();
  associatedlc_to_cp.clear();
  associatedcp_to_simlc.clear();
  associatedsimlc_to_cp.clear();
  caloparticle_hit_eta.clear();
  caloparticle_hit_phi.clear();
  caloparticle_hit_energy.clear();
  caloparticle_hit_event.clear();
  caloparticle_event.clear();
  shared_energy.clear();
  layercluster_hit_eta.clear();
  layercluster_hit_phi.clear();
  layercluster_hit_density.clear();
  layercluster_hit_energy.clear();
  layercluster_hit_event.clear();
  layercluster_event.clear();
  layercluster_nhits.clear();
  layercluster_resolution.clear();
  pfrechit_eta.clear();
  pfrechit_phi.clear();
  pfrechit_energy.clear();
  pfrechit_time.clear();
  simcl_energy.clear();
  simcl_eta.clear();
  simcl_phi.clear();
  sc2lc_score.clear();
  simcl_layer.clear();
  associatedsc_to_lc.clear();
  lc2sc_score.clear();
  associatedlc_to_sc.clear();
  simcl_shared_energy.clear();
  caloparticle_bx.clear();

  seed_density.clear();
  simlayercluster_energy.clear();
  simlayercluster_eta.clear();
  simlayercluster_phi.clear();
  simlayercluster_seed_eta.clear();
  simlayercluster_seed_phi.clear();
  simlayercluster_hit_eta.clear();
  simlayercluster_hit_phi.clear();
  simlayercluster_hit_energy.clear();
  simlayercluster_nhits.clear();
  simlayercluster_resolution.clear();
}

float LayerClusterDumper::distance (GlobalPoint point1, GlobalPoint point2) {
  float deltaEta = point1.eta() - point2.eta();
  auto o2pi = 1./(2*M_PI);
  float deltaPhi = point1.phi() - point2.phi();
  if (std::abs(deltaPhi) > M_PI) {
    auto n = std::round(deltaPhi*o2pi);
    deltaPhi = deltaPhi - n*M_PI;
  }
  return std::sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

std::vector<std::pair<std::pair<float, float>, float>> LayerClusterDumper::addHitsAndDensities(std::vector<std::pair<uint32_t, float>> hits_and_energies, const CaloGeometry* geom, float totalEnergy) {
  std::vector<std::pair<std::pair<float, float>, float>> hits_and_densities;
  hits_and_densities.resize(hits_and_energies.size());
  for (const auto& hit_and_energy : hits_and_energies) {
    uint32_t hit = hit_and_energy.first;
    float energy = hit_and_energy.second;
    GlobalPoint position = geom->getPosition(DetId(hit));
    float density = energy;
    for (const auto& other_hit_and_energy : hits_and_energies) {
      uint32_t other_hit = other_hit_and_energy.first;
      if (other_hit == hit ) continue;
      float other_energy = other_hit_and_energy.second;
      GlobalPoint other_position = geom->getPosition(DetId(other_hit));
      if (distance(position, other_position) < 1.5*0.0175)
	density += 0.5*other_energy;
    }
    //density *= totalEnergy;
    hits_and_densities.push_back(std::make_pair(std::make_pair(position.eta(), position.phi()), density));
 }
 return hits_and_densities;
}

std::vector<std::pair<std::pair<float, float>, float>> LayerClusterDumper::addHitsAndDensities(std::vector<std::pair<DetId, float>> hits_and_fractions, const CaloGeometry* geom, float totalEnergy) {
  std::vector<std::pair<std::pair<float, float>, float>> hits_and_densities;
  hits_and_densities.resize(hits_and_fractions.size());
  for (const auto& hit_and_fraction : hits_and_fractions) {
    DetId hit = hit_and_fraction.first;
    float fraction = hit_and_fraction.second;
    GlobalPoint position = geom->getPosition(hit);
    float density = fraction;
    for (const auto& other_hit_and_fraction : hits_and_fractions) {
      DetId other_hit = other_hit_and_fraction.first;
      if (other_hit == hit ) continue;
      float other_fraction = other_hit_and_fraction.second;
      GlobalPoint other_position = geom->getPosition(other_hit);
      if (distance(position, other_position) < 1.5*0.0175)
	density += 0.5*other_fraction;
    }
    density *= totalEnergy;
    hits_and_densities.push_back(std::make_pair(std::make_pair(position.eta(), position.phi()), density));
 }
 return hits_and_densities;
}

LayerClusterDumper::LayerClusterDumper(const edm::ParameterSet& ps)
  : pfrechit_token_(consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("pfrechits"))),
    layer_clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layerclusters"))),
    sim_layer_clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("simlayerclusters"))),
    caloparticles_token_(consumes<std::vector<CaloParticle>>(ps.getParameter<edm::InputTag>("caloparticles"))), 
    simclusters_token_(consumes<std::vector<SimCluster>>(ps.getParameter<edm::InputTag>("caloparticles"))), 
    simtoreco_token_(consumes<hgcal::SimToRecoCollection>(ps.getParameter<edm::InputTag>("simToRecoCollection"))),
    recotosim_token_(consumes<hgcal::RecoToSimCollection>(ps.getParameter<edm::InputTag>("recoToSimCollection"))),
    simtorecosim_token_(consumes<hgcal::SimToRecoCollection>(ps.getParameter<edm::InputTag>("simToRecoCollectionSim"))),
    recotosimsim_token_(consumes<hgcal::RecoToSimCollection>(ps.getParameter<edm::InputTag>("recoToSimCollectionSim"))),
    simtoreco_simcl_token_(consumes<hgcal::SimToRecoCollectionWithSimClusters>(ps.getParameter<edm::InputTag>("simToRecoCollectionSC"))),
    recotosim_simcl_token_(consumes<hgcal::RecoToSimCollectionWithSimClusters>(ps.getParameter<edm::InputTag>("recoToSimCollectionSC"))),
    geometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()) {
};

LayerClusterDumper::~LayerClusterDumper() { clearVariables(); };

void LayerClusterDumper::beginJob() {
  edm::Service<TFileService> fs;
  pfrechit_tree_ = fs->make<TTree>("pfrechits", "PFRecHits");
  layercluster_tree_ = fs->make<TTree>("layerclusters", "Layer Clusters");
  simlayercluster_tree_ = fs->make<TTree>("simlayerclusters", "SimLayerClusters");
  caloparticle_tree_ = fs->make<TTree>("caloparticles", "CaloParticles");
  simcluster_tree_ = fs->make<TTree>("simclusters", "SimClusters");

  pfrechit_tree_->Branch("pfrechitEta", &pfrechit_eta);
  pfrechit_tree_->Branch("pfrechitPhi", &pfrechit_phi);
  pfrechit_tree_->Branch("pfrechitEnergy", &pfrechit_energy);
  pfrechit_tree_->Branch("pfrechitTime", &pfrechit_time);
  layercluster_tree_->Branch("layerClusterEnergy", &layercluster_energy);
  layercluster_tree_->Branch("layerClusterEta", &layercluster_eta);
  layercluster_tree_->Branch("layerClusterPhi", &layercluster_phi);
  layercluster_tree_->Branch("layerClusterHitEta", &layercluster_hit_eta);
  layercluster_tree_->Branch("layerClusterHitPhi", &layercluster_hit_phi);
  layercluster_tree_->Branch("layerClusterHitDensity", &layercluster_hit_density);
  layercluster_tree_->Branch("layerClusterHitEnergy", &layercluster_hit_energy);
  layercluster_tree_->Branch("layerClusterEvent", &layercluster_event);
  layercluster_tree_->Branch("layerClusterNumberOfHits", &layercluster_nhits);
  layercluster_tree_->Branch("seedDensity", &seed_density);
  layercluster_tree_->Branch("layerClusterResolution", &layercluster_resolution);
  layercluster_tree_->Branch("layerClusterSeedEta", &layercluster_seed_eta);
  layercluster_tree_->Branch("layerClusterSeedPhi", &layercluster_seed_phi);
  layercluster_tree_->Branch("layerClusterLayer", &layercluster_layer);
  layercluster_tree_->Branch("layerClusterPUContribution", &layercluster_pu_contribution);
  simlayercluster_tree_->Branch("simLayerClusterEnergy", &simlayercluster_energy);
  simlayercluster_tree_->Branch("simLayerClusterEta", &simlayercluster_eta);
  simlayercluster_tree_->Branch("simLayerClusterPhi", &simlayercluster_phi);
  simlayercluster_tree_->Branch("simLayerClusterHitEnergy", &simlayercluster_hit_energy);
  simlayercluster_tree_->Branch("simLayerClusterHitEta", &simlayercluster_hit_eta);
  simlayercluster_tree_->Branch("simLayerClusterHitPhi", &simlayercluster_hit_phi);
  simlayercluster_tree_->Branch("simLayerClusterResolution", &simlayercluster_resolution);
  simlayercluster_tree_->Branch("simLayerClusterSeedEta", &simlayercluster_seed_eta);
  simlayercluster_tree_->Branch("simLayerClusterSeedPhi", &simlayercluster_seed_phi);

  caloparticle_tree_->Branch("caloParticleEnergy", &caloparticle_energy);
  caloparticle_tree_->Branch("caloParticleEta", &caloparticle_eta);
  caloparticle_tree_->Branch("caloParticlePhi", &caloparticle_phi);
  caloparticle_tree_->Branch("caloParticleHitEta", &caloparticle_hit_eta);
  caloparticle_tree_->Branch("caloParticleHitPhi", &caloparticle_hit_phi);
  caloparticle_tree_->Branch("caloParticleHitDensity", &caloparticle_hit_density);
  caloparticle_tree_->Branch("caloParticleHitEnergy", &caloparticle_hit_energy);
  caloparticle_tree_->Branch("caloParticleEvent", &caloparticle_event);
  caloparticle_tree_->Branch("caloParticleBX", &caloparticle_bx); 
  caloparticle_tree_->Branch("sharedEnergy", &shared_energy);

  simcluster_tree_->Branch("simClusterEnergy", &simcl_energy);
  simcluster_tree_->Branch("simClusterEta", &simcl_eta);
  simcluster_tree_->Branch("simClusterPhi", &simcl_phi);
  simcluster_tree_->Branch("simClusterLayer", &simcl_layer);
  simcluster_tree_->Branch("sharedEnergy", &simcl_shared_energy);
  simcluster_tree_->Branch("simToRecoAssociation", &sc2lc_score);
  layercluster_tree_->Branch("recoToSimAssociation", &lc2sc_score);
  simcluster_tree_->Branch("AssociatedLC", &associatedlc_to_sc);
  layercluster_tree_->Branch("AssociatedSC", &associatedsc_to_lc);

  caloparticle_tree_->Branch("simToRecoAssociation", &cp2lc_score);
  caloparticle_tree_->Branch("recoToSimAssociationSim", &cp2lc_score_sim);
  caloparticle_tree_->Branch("AssociatedLC", &associatedlc_to_cp);
  caloparticle_tree_->Branch("AssociatedSimLC", &associatedsimlc_to_cp);
  layercluster_tree_->Branch("recoToSimAssociation", &lc2cp_score);
  layercluster_tree_->Branch("AssociatedCP", &associatedcp_to_lc);
  simlayercluster_tree_->Branch("recoToSimAssociationSim", &lc2cp_score_sim);
  simlayercluster_tree_->Branch("AssociatedCP", &associatedcp_to_simlc);
  event_index = 0;
}

void LayerClusterDumper::beginRun(edm::Run const&, edm::EventSetup const& es) {
  geom = &(es.getData(geometry_token_));
}

void LayerClusterDumper::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  event_index++;
  clearVariables();

  edm::Handle<reco::PFRecHitCollection> pfrechits_h;
  event.getByToken(pfrechit_token_, pfrechits_h);
  const auto& pfrechits = *pfrechits_h;
  
  edm::Handle<std::vector<reco::CaloCluster>> layer_clusters_h;
  event.getByToken(layer_clusters_token_, layer_clusters_h);
  const auto& layer_clusters = *layer_clusters_h;

  edm::Handle<std::vector<reco::CaloCluster>> sim_layer_clusters_h;
  event.getByToken(sim_layer_clusters_token_, sim_layer_clusters_h);
  const auto& sim_layer_clusters = *sim_layer_clusters_h;

  edm::Handle<std::vector<CaloParticle>> caloparticle_h;
  event.getByToken(caloparticles_token_, caloparticle_h);
  const auto& caloparticles = *caloparticle_h;

  edm::Handle<std::vector<SimCluster>> simcluster_h;
  event.getByToken(simclusters_token_, simcluster_h);
  const auto& simclusters = *simcluster_h;

  edm::Handle<hgcal::SimToRecoCollection> simToReco_h;
  event.getByToken(simtoreco_token_, simToReco_h);
  const auto& simToReco = *simToReco_h;

  edm::Handle<hgcal::RecoToSimCollection> recoToSim_h;
  event.getByToken(recotosim_token_, recoToSim_h);
  const auto& recoToSim = *recoToSim_h;
 
  edm::Handle<hgcal::SimToRecoCollection> simToRecoSim_h;
  event.getByToken(simtorecosim_token_, simToRecoSim_h);
  const auto& simToRecoSim = *simToRecoSim_h;

  edm::Handle<hgcal::RecoToSimCollection> recoToSimSim_h;
  event.getByToken(recotosimsim_token_, recoToSimSim_h);
  const auto& recoToSimSim = *recoToSimSim_h;

  edm::Handle<hgcal::SimToRecoCollectionWithSimClusters> simToReco_scl_h;
  event.getByToken(simtoreco_simcl_token_, simToReco_scl_h);
  const auto& simToReco_simcl = *simToReco_scl_h;
						       
  edm::Handle<hgcal::RecoToSimCollectionWithSimClusters> recoToSim_scl_h;
  event.getByToken(recotosim_simcl_token_, recoToSim_scl_h);
  const auto& recoToSim_simcl = *recoToSim_scl_h;

  ev_event_ = event_index;
  nclusters_ = layer_clusters.size();
  nsimlayerclusters_ = sim_layer_clusters.size();

  ncaloparticles_ = caloparticles.size(); 
  npfrechits_ = pfrechits.size();
  nsimclusters_ = simclusters.size();
  /*pfrechit_eta.resize(npfrechits_);
  pfrechit_phi.resize(npfrechits_);
  pfrechit_energy.resize(npfrechits_);
  pfrechit_time.resize(npfrechits_);*/

  layercluster_hit_energy.resize(nclusters_);
  layercluster_hit_eta.resize(nclusters_);
  layercluster_hit_phi.resize(nclusters_);
  layercluster_hit_density.resize(nclusters_);
  layercluster_resolution.resize(nclusters_);
  layercluster_layer.resize(nclusters_);
  seed_density.resize(nclusters_);
  layercluster_seed_eta.resize(nclusters_);
  layercluster_seed_phi.resize(nclusters_);

  simlayercluster_hit_energy.resize(nsimlayerclusters_);
  simlayercluster_hit_eta.resize(nsimlayerclusters_);
  simlayercluster_hit_phi.resize(nsimlayerclusters_);
  simlayercluster_resolution.resize(nsimlayerclusters_);
  simlayercluster_seed_eta.resize(nsimlayerclusters_);
  simlayercluster_seed_phi.resize(nsimlayerclusters_);

  caloparticle_hit_energy.resize(ncaloparticles_);
  caloparticle_hit_eta.resize(ncaloparticles_);
  caloparticle_hit_phi.resize(ncaloparticles_);
  caloparticle_hit_density.resize(ncaloparticles_);
  shared_energy.resize(ncaloparticles_);

  std::unordered_map<DetId, float> hitmap;
  for (auto pfrechits_iterator = pfrechits.begin(); pfrechits_iterator != pfrechits.end(); ++pfrechits_iterator) {
    DetId id(pfrechits_iterator->detId());
    //GlobalPoint position = geom->getPosition(id);
    //pfrechit_eta.push_back(position.eta());
    //pfrechit_phi.push_back(position.phi());
    float energy = pfrechits_iterator->energy();
    //pfrechit_energy.push_back(energy);
    //pfrechit_time.push_back(pfrechits_iterator->time());
    hitmap.insert(std::make_pair(id, energy));
  }

  std::vector<size_t> cPIndices;
  removeCPFromPU(caloparticles, cPIndices);

  std::unordered_map<DetId, std::vector<std::pair<size_t, float>>> detId2cpMap;
  size_t cpIndex = 0;
  for (auto cp_iterator = caloparticles.begin(); cp_iterator != caloparticles.end(); ++cp_iterator) {
    std::vector<std::pair<uint32_t, float>> haf = cp_iterator->hits_and_fractions();
    for (const auto& cl : cp_iterator->simClusters()) {
      const std::vector<std::pair<uint32_t, float>> haf = cl->hits_and_fractions();
      for (const auto& hit : haf) {
	DetId id(hit.first);
	float fraction = hit.second;
	if (detId2cpMap.find(id) == detId2cpMap.end()) {
	  std::vector<std::pair<size_t, float>> idxAndFraction{std::make_pair(cpIndex, fraction)};
	  detId2cpMap.insert(std::make_pair(id, idxAndFraction));
	} else {
	    detId2cpMap[id].push_back(std::make_pair(cpIndex, fraction));
	}
      }
    }
    cpIndex++;
  }

  lc2cp_score.resize(nclusters_);
  associatedcp_to_lc.resize(nclusters_);
  lc2sc_score.resize(nclusters_);
  associatedsc_to_lc.resize(nclusters_);
  int lc_index = 0;
  for (auto lc_iterator = layer_clusters.begin(); lc_iterator != layer_clusters.end(); ++lc_iterator) {
    layercluster_energy.push_back(lc_iterator->energy());
    layercluster_eta.push_back(lc_iterator->eta());
    layercluster_phi.push_back(lc_iterator->phi());
    DetId seedId = lc_iterator->seed();
    int layerId = 0;
    if (seedId.det() == DetId::Hcal) {
      HcalDetId hid(seedId);
      layerId = hid.depth();
      if (seedId.subdetId() == HcalSubdetector::HcalOuter)
	layerId += 1;
    }
    layercluster_layer.push_back(layerId);
    float seedEta = (geom->getPosition(seedId)).eta();
    float seedPhi = (geom->getPosition(seedId)).phi();
    layercluster_seed_eta.push_back(seedEta);
    layercluster_seed_phi.push_back(seedPhi);
    std::vector<std::pair<DetId, float>> hits_and_fractions = lc_iterator->hitsAndFractions();
    layercluster_nhits.push_back(hits_and_fractions.size());
    float pu_contribution_to_cluster = 0.f;
    for (const auto& haf : hits_and_fractions) { 
      float pu_contribution_to_hit = 0.f;
      //float eta = (geom->getPosition(haf.first)).eta();
      //float phi = (geom->getPosition(haf.second)).phi();
      //layercluster_hit_eta[lc_index].push_back(eta);
      //layercluster_hit_phi[lc_index].push_back(phi);
      layercluster_hit_energy[lc_index].push_back(haf.second*lc_iterator->energy());
      if (detId2cpMap.find(haf.first) == detId2cpMap.end()) continue;
      auto cp_and_fraction = detId2cpMap[haf.first];
      for (const auto& [cp, frac] : cp_and_fraction) {
	if (std::find(cPIndices.begin(), cPIndices.end(), cp) != cPIndices.end()) continue;
	pu_contribution_to_hit += frac;
      }
      float hit_energy = hitmap[haf.first];
      pu_contribution_to_hit *= hit_energy;
      pu_contribution_to_cluster += pu_contribution_to_hit;
    }
    float lc_energy = lc_iterator->energy();
    layercluster_pu_contribution.push_back(pu_contribution_to_cluster/lc_energy);
    std::vector<std::pair<std::pair<float, float>, float>> hits_and_densities = addHitsAndDensities(hits_and_fractions, geom, lc_energy);
    for (const auto& hit_and_density : hits_and_densities) {
      //DetId hit = hit_and_density.first;
      float eta = hit_and_density.first.first;
      float phi = hit_and_density.first.second;
      float density = hit_and_density.second;
      if (eta == seedEta && phi == seedPhi) 
	seed_density[lc_index].push_back(density);
//      GlobalPoint position = geom->getPosition(hit);
      layercluster_hit_eta[lc_index].push_back(eta);
      layercluster_hit_phi[lc_index].push_back(phi);
      layercluster_hit_density[lc_index].push_back(density);
      layercluster_hit_event.push_back(event_index);
    }
    layercluster_event.push_back(event_index);
    const edm::Ref<std::vector<reco::CaloCluster>> lcRef(layer_clusters_h, lc_index);
    const auto& cpsIt = recoToSim.find(lcRef);
    if (cpsIt != recoToSim.end()) {
      const auto& cps = cpsIt->val;
      for (const auto& cpPair : cps) {
	auto cp_id = (cpPair.first).get() - (edm::Ref<std::vector<CaloParticle>>(caloparticle_h, 0)).get();
	associatedcp_to_lc[lc_index].push_back(cp_id);
	lc2cp_score[lc_index].push_back(cpPair.second);
      }
    }
    const auto& scsIt = recoToSim_simcl.find(lcRef);
    if (scsIt != recoToSim_simcl.end()) {
      const auto& scs = scsIt->val;
      for (const auto& scPair : scs) {
	auto sc_id = (scPair.first).get() - (edm::Ref<std::vector<SimCluster>>(simcluster_h, 0)).get();
	associatedsc_to_lc[lc_index].push_back(sc_id);
	lc2sc_score[lc_index].push_back(scPair.second);
      }
    }
    lc_index++;
  }


  cp2lc_score.resize(caloparticles.size());
  cp2lc_score_sim.resize(caloparticles.size());
  associatedlc_to_cp.resize(caloparticles.size());
  associatedsimlc_to_cp.resize(caloparticles.size());
  int cp_index = 0;
  for (auto cp_iterator = caloparticles.begin(); cp_iterator != caloparticles.end(); ++cp_iterator) {
    float cp_energy = cp_iterator->energy();
    caloparticle_energy.push_back(cp_energy);
    caloparticle_eta.push_back(cp_iterator->eta());
    caloparticle_phi.push_back(cp_iterator->phi());
    caloparticle_bx.push_back(cp_iterator->g4Tracks()[0].eventId().bunchCrossing());
    std::vector<std::pair<std::pair<float, float>, float>> hits_and_energies;
    for (const auto& cl : cp_iterator->simClusters()) {
      const std::vector<std::pair<uint32_t, float>> hae = (*cl).hits_and_energies();
      for (auto it = hae.begin(); it != hae.end(); ++it) {
	//if (it->second < 0.1) continue;
	//caloparticle_hit_energy.push_back(it->second);
	DetId id(it->first);
	float eta = (geom->getPosition(id)).eta();
	float phi = (geom->getPosition(id)).phi();
	hits_and_energies.push_back(std::make_pair(std::make_pair(eta,phi), it->second));
      }
    }
    //std::vector<std::pair<std::pair<float, float>, float>> hits_and_densities = addHitsAndDensities(hits_and_energies, geom, cp_energy);
    //for (const auto& hit_and_density : hits_and_densities) {
      //float eta = hit_and_density.first.first;
      //float phi = hit_and_density.first.second;
      //float density = hit_and_density.second;
      //caloparticle_hit_eta.push_back(eta);
      //caloparticle_hit_phi.push_back(phi);
      //caloparticle_hit_density.push_back(density);
      //caloparticle_hit_event.push_back(event_index);
    //}
    //for (const auto& hit_and_energy : hits_and_energies) {
    //  caloparticle_hit_energy[cp_index].push_back(hit_and_energy.second);
    //  caloparticle_hit_eta[cp_index].push_back(hit_and_energy.first.first);
    //  caloparticle_hit_phi[cp_index].push_back(hit_and_energy.first.second);
    //}
    caloparticle_event.push_back(cp_iterator->g4Tracks()[0].eventId().event());
    const edm::Ref<std::vector<CaloParticle>> cpRef(caloparticle_h, cp_index);
    const auto& lcsIt = simToReco.find(cpRef);
    const auto& simlcsIt = simToRecoSim.find(cpRef);
    if (lcsIt != simToReco.end()) {
      const auto& lcs = lcsIt->val;
      for (const auto& lcPair : lcs) {
	auto lc_id = (lcPair.first).get() - (edm::Ref<std::vector<reco::CaloCluster>>(layer_clusters_h, 0)).get();
	associatedlc_to_cp[cp_index].push_back(lc_id);
	cp2lc_score[cp_index].push_back(lcPair.second.second);
	shared_energy[cp_index].push_back(lcPair.second.first);
      }
    }
    if (simlcsIt != simToRecoSim.end()) {
      const auto& simlcs = simlcsIt->val;
      for (const auto& simlcPair : simlcs) {
	auto simlc_id = (simlcPair.first).get() - (edm::Ref<std::vector<reco::CaloCluster>>(sim_layer_clusters_h, 0)).get();
	associatedsimlc_to_cp[cp_index].push_back(simlc_id);
	cp2lc_score_sim[cp_index].push_back(simlcPair.second.second);
      }
    }
    cp_index++;
  }

  sc2lc_score.resize(nsimclusters_);
  associatedlc_to_sc.resize(nsimclusters_);
  auto sc_index = 0;
  simcl_shared_energy.resize(nsimclusters_);
  for (auto sc_iterator = simclusters.begin(); sc_iterator != simclusters.end(); ++sc_iterator) {
    simcl_energy.push_back(sc_iterator->energy());
    simcl_eta.push_back(sc_iterator->eta());
    simcl_phi.push_back(sc_iterator->phi());
    auto firstHitId = ((sc_iterator->hits_and_fractions())[0]).first;
    HcalDetId id(firstHitId);
    auto layer = (id.subdetId() == HcalBarrel) ? id.depth() : id.depth() + 1;
    simcl_layer.push_back(layer);

    const edm::Ref<std::vector<SimCluster>> scRef(simcluster_h, sc_index);
    const auto& lcsIt = simToReco_simcl.find(scRef);
    if (lcsIt != simToReco_simcl.end())  {
      const auto& lcs = lcsIt->val;
      for (const auto& lcPair : lcs) {
	auto lc_id = (lcPair.first).get() - (edm::Ref<std::vector<reco::CaloCluster>>(layer_clusters_h, 0)).get();
	associatedlc_to_sc[sc_index].push_back(lc_id);
	sc2lc_score[sc_index].push_back(lcPair.second.second);
	simcl_shared_energy[sc_index].push_back(lcPair.second.first);
      }
    }
    sc_index++;
  }

  lc2cp_score_sim.resize(nsimlayerclusters_);
  associatedcp_to_simlc.resize(nsimlayerclusters_);
  int simlc_index = 0;
  for (auto simlc_iterator = sim_layer_clusters.begin(); simlc_iterator != sim_layer_clusters.end(); ++simlc_iterator) {
    simlayercluster_energy.push_back(simlc_iterator->energy());
    simlayercluster_phi.push_back(simlc_iterator->phi());
    simlayercluster_eta.push_back(simlc_iterator->eta());
    DetId seedId = simlc_iterator->seed();
    float seedEta = (geom->getPosition(seedId)).eta();
    float seedPhi = (geom->getPosition(seedId)).phi();
    simlayercluster_seed_eta.push_back(seedEta);
    simlayercluster_seed_phi.push_back(seedPhi);
    std::vector<std::pair<DetId, float>> hitsAndFractions = simlc_iterator->hitsAndFractions();
    //for (const auto& haf : hitsAndFractions) {
      //GlobalPoint pos = geom->getPosition(haf.first);
      //float eta = pos.eta();
      //float phi = pos.phi();
      //simlayercluster_hit_eta[simlc_index].push_back(eta);
      //simlayercluster_hit_phi[simlc_index].push_back(phi);
      //simlayercluster_hit_energy[simlc_index].push_back(haf.second * simlc_iterator->energy());
    //}
    const edm::Ref<std::vector<reco::CaloCluster>> simlcRef(sim_layer_clusters_h, simlc_index);
    const auto& cpsIt = recoToSimSim.find(simlcRef);
    if (cpsIt != recoToSimSim.end()) {
      const auto& cps = cpsIt->val;
      for (const auto& cpPair : cps) {
	auto cp_id = (cpPair.first).get() - (edm::Ref<std::vector<CaloParticle>>(caloparticle_h, 0)).get();
	associatedcp_to_simlc[simlc_index].push_back(cp_id);
	lc2cp_score_sim[simlc_index].push_back(cpPair.second);
      }
    }                                                                                                                    
    simlc_index++;
  }
  pfrechit_tree_->Fill();
  layercluster_tree_->Fill();
  simlayercluster_tree_->Fill();
  simcluster_tree_->Fill();
  caloparticle_tree_->Fill();
}

void LayerClusterDumper::endJob() {}

void LayerClusterDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pfrechits", edm::InputTag("particleFlowRecHitECAL", ""));
  desc.add<edm::InputTag>("layerclusters", edm::InputTag("barrelLayerClusters"));
  desc.add<edm::InputTag>("simlayerclusters", edm::InputTag("simBarrelLayerClusters"));
  desc.add<edm::InputTag>("caloparticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("simclusters", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("simToRecoCollection", edm::InputTag("barrelLayerClusterCaloParticleAssociationProducer"));
  desc.add<edm::InputTag>("recoToSimCollection", edm::InputTag("barrelLayerClusterCaloParticleAssociationProducer"));
  desc.add<edm::InputTag>("simToRecoCollectionSim", edm::InputTag("simBarrelLayerClusterCaloParticleAssociationProducer"));
  desc.add<edm::InputTag>("recoToSimCollectionSim", edm::InputTag("simBarrelLayerClusterCaloParticleAssociationProducer"));
  desc.add<edm::InputTag>("simToRecoCollectionSC", edm::InputTag("barrelLayerClusterSimClusterAssociationProducer"));
  desc.add<edm::InputTag>("recoToSimCollectionSC", edm::InputTag("barrelLayerClusterSimClusterAssociationProducer"));
  descriptions.add("layerClusterDumper", desc);
}

DEFINE_FWK_MODULE(LayerClusterDumper);
  
