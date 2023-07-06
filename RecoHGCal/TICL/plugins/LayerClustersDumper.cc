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
#include "DataFormats/Math/interface/deltaPhi.h"
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
#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociatorBaseImpl.h"
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
  typedef std::vector<std::pair<uint32_t, float>> hits_and_floats;
  typedef std::vector<std::pair<DetId, float>> ids_and_floats;

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override {};

  //void initialize();
  //void buildLayers();

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override;

  //LCs
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_token_;
  //CPs
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;
  //SCs
  //const edm::EDGetTokenT<std::vector<SimCluster>> simclusters_token_;
  //Association maps
  const edm::EDGetTokenT<hgcal::SimToRecoCollection> simtoreco_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollection> recotosim_token_;

  //Geometry
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  
  void clearVariables();

  template<typename T, typename V>
  void addHitsAndDensities(V& hits_and_densities, T& hits_and_fractions, const CaloGeometry* geom, float totalEnergy);
  
  float distance(GlobalPoint point1, GlobalPoint point2);
  unsigned int event_index;
  TTree* tree_;

  unsigned int ev_event_;
  unsigned int nclusters_;

  std::vector<float> layercluster_energy;
  std::vector<float> layercluster_eta;
  std::vector<float> layercluster_phi;
  std::vector<float> layercluster_hit_eta;
  std::vector<float> layercluster_hit_phi;
  std::vector<float> layercluster_hit_density;

  std::vector<float> caloparticle_energy;
  std::vector<float> caloparticle_eta;
  std::vector<float> caloparticle_phi;
  std::vector<float> caloparticle_hit_eta;
  std::vector<float> caloparticle_hit_phi;
  std::vector<float> caloparticle_hit_density;

  std::vector<std::vector<float>> cp2lc_score;
  std::vector<std::vector<float>> lc2cp_score;

  TTree* layercluster_tree_;
  TTree* caloparticle_tree_;
  TTree* association_tree_;

};

void LayerClusterDumper::clearVariables() {
  layercluster_energy.clear();
  layercluster_eta.clear();
  layercluster_phi.clear();
  caloparticle_energy.clear();
  caloparticle_eta.clear();
  caloparticle_phi.clear();
  cp2lc_score.clear();
  lc2cp_score.clear();
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

template<typename T, typename V>
void LayerClusterDumper::addHitsAndDensities(V& hits_and_densities, T& hits_and_fractions, const CaloGeometry* geom, float totalEnergy) {
  for (const auto& hit_and_fraction : hits_and_fractions) {
    auto hit = hit_and_fraction.first;
    float fraction = hit_and_fraction.second;
    GlobalPoint position = geom->getPosition(hit);
    float density = fraction;
    for (const auto& other_hit_and_fraction : hits_and_fractions) {
      auto other_hit = other_hit_and_fraction.first;
      if (other_hit == hit ) continue;
      float other_fraction = other_hit_and_fraction.second;
      GlobalPoint other_position = geom->getPosition(other_hit);
      if (distance(position, other_position) < 1.5*0.0175)
	density += 0.5*other_fraction;
    }
    density *= totalEnergy;
    hits_and_densities.push_back(std::make_pair(hit, density));
 }
}

LayerClusterDumper::LayerClusterDumper(const edm::ParameterSet& ps)
  : layer_clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layerclusters"))),
    caloparticles_token_(consumes<std::vector<CaloParticle>>(ps.getParameter<edm::InputTag>("caloparticles"))),
    simtoreco_token_(consumes<hgcal::SimToRecoCollection>(ps.getParameter<edm::InputTag>("simToRecoCollection"))),
    recotosim_token_(consumes<hgcal::RecoToSimCollection>(ps.getParameter<edm::InputTag>("recoToSimCollection"))) {
      
      geometry_token_ = esConsumes<CaloGeometry, CaloGeometryRecord>();
};

LayerClusterDumper::~LayerClusterDumper() { clearVariables(); };

void LayerClusterDumper::beginJob() {
  edm::Service<TFileService> fs;
  layercluster_tree_ = fs->make<TTree>("layerclusters", "Layer Clusters");
  association_tree_ = fs->make<TTree>("associartions", "Associations");
  caloparticle_tree_ = fs->make<TTree>("caloparticles", "CaloParticles");

  layercluster_tree_->Branch("layerClusterEnergy", &layercluster_energy);
  layercluster_tree_->Branch("layerClusterEta", &layercluster_eta);
  layercluster_tree_->Branch("layerClusterPhi", &layercluster_phi);
  layercluster_tree_->Branch("layerClusterHitEta", &layercluster_hit_eta);
  layercluster_tree_->Branch("layerClusterHitPhi", &layercluster_hit_phi);
  layercluster_tree_->Branch("layerClusterHitDensity", &layercluster_hit_density);

  caloparticle_tree_->Branch("caloParticleEnergy", &caloparticle_energy);
  caloparticle_tree_->Branch("caloParticleEta", &caloparticle_eta);
  caloparticle_tree_->Branch("caloParticlePhi", &caloparticle_phi);
  caloparticle_tree_->Branch("caloParticleHitEta", &caloparticle_hit_eta);
  caloparticle_tree_->Branch("caloParticleHitPhi", &caloparticle_hit_phi);
  caloparticle_tree_->Branch("caloParticleHitDensity", &caloparticle_hit_density);

  association_tree_->Branch("simToRecoAssociation", &cp2lc_score);
  association_tree_->Branch("recoToSimAssociation", &lc2cp_score);

  event_index = 0;
}

void LayerClusterDumper::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  event_index++;
  clearVariables();

  const CaloGeometry* geom = &(setup.getData(geometry_token_));

  edm::Handle<std::vector<reco::CaloCluster>> layer_clusters_h;
  event.getByToken(layer_clusters_token_, layer_clusters_h);
  const auto& layer_clusters = *layer_clusters_h;

  edm::Handle<std::vector<CaloParticle>> caloparticle_h;
  event.getByToken(caloparticles_token_, caloparticle_h);
  const auto& caloparticles = *caloparticle_h;

  edm::Handle<hgcal::SimToRecoCollection> simToReco_h;
  event.getByToken(simtoreco_token_, simToReco_h);
  const auto& simToReco = *simToReco_h;

  edm::Handle<hgcal::RecoToSimCollection> recoToSim_h;
  event.getByToken(recotosim_token_, recoToSim_h);
  const auto& recoToSim = *recoToSim_h;
  
  ev_event_ = event_index;
  nclusters_ = layer_clusters.size();
  

  for (auto lc_iterator = layer_clusters.begin(); lc_iterator != layer_clusters.end(); ++lc_iterator) {
    layercluster_energy.push_back(lc_iterator->energy());
    layercluster_eta.push_back(lc_iterator->eta());
    layercluster_phi.push_back(lc_iterator->phi());
    ids_and_floats hits_and_fractions = lc_iterator->hitsAndFractions();
    float lc_energy = lc_iterator->energy();
    ids_and_floats hits_and_densities;
    hits_and_densities.resize(hits_and_fractions.size());
    addHitsAndDensities(hits_and_densities, hits_and_fractions, geom, lc_energy);
    for (const auto& hit_and_density : hits_and_densities) {
      uint32_t hit = hit_and_density.first;
      float density = hit_and_density.second;
      GlobalPoint position = geom->getPosition(hit);
      layercluster_hit_eta.push_back(position.eta());
      layercluster_hit_phi.push_back(position.phi());
      layercluster_hit_density.push_back(density);
    }
  }

  for (auto cp_iterator = caloparticles.begin(); cp_iterator != caloparticles.end(); ++cp_iterator) {
    caloparticle_energy.push_back(cp_iterator->energy());
    caloparticle_eta.push_back(cp_iterator->eta());
    caloparticle_phi.push_back(cp_iterator->phi());
    const hits_and_floats hits_and_fractions = cp_iterator->hits_and_fractions();
    float cp_energy = cp_iterator->energy();
    hits_and_floats hits_and_densities;
    hits_and_densities.resize(hits_and_fractions.size());
    addHitsAndDensities(hits_and_densities,hits_and_fractions, geom, cp_energy);
    for (const auto& hit_and_density : hits_and_densities) {
      uint32_t hit = hit_and_density.first;
      float density = hit_and_density.second;
      GlobalPoint position = geom->getPosition(hit);
      caloparticle_hit_eta.push_back(position.eta());
      caloparticle_hit_phi.push_back(position.phi());
      caloparticle_hit_density.push_back(density);
    }
}

  lc2cp_score.resize(layer_clusters.size());
  for (unsigned int lcId = 0; lcId < layer_clusters.size(); ++lcId) {
    const edm::Ref<std::vector<reco::CaloCluster>> lcRef(layer_clusters_h, lcId);
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
    const auto& lcs = lcsIt->val;
    for (const auto& lcPair : lcs) {
      cp2lc_score[cpId].push_back(lcPair.second.second);
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
  desc.add<edm::InputTag>("simToRecoCollection", edm::InputTag("barrelLayerClusterCaloParticleAssociationProducer"));
  desc.add<edm::InputTag>("recoToSimCollection", edm::InputTag("barrelLayerClusterCaloParticleAssociationProducer"));

  descriptions.add("layerClusterDumper", desc);
}

DEFINE_FWK_MODULE(LayerClusterDumper);
  
