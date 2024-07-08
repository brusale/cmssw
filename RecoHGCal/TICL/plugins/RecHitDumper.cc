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
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
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

class RecHitDumper : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
public:
  explicit RecHitDumper(const edm::ParameterSet&);
  ~RecHitDumper() override;
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
  const edm::EDGetTokenT<EcalRecHitCollection> ecalrechit_token_;
  const edm::EDGetTokenT<edm::PCaloHitContainer> simhit_token_;
  //Geometry
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  
  void clearVariables();


  TTree* tree_;
  
  const CaloGeometry* geom;

  unsigned int ev_event_;
  unsigned int npfrechits_;

  std::vector<float> pfrechit_eta;
  std::vector<float> pfrechit_phi;
  std::vector<float> pfrechit_energy;
  std::vector<float> pfrechit_time;
  std::vector<int> pfrechit_ieta;
  std::vector<int> pfrechit_iphi;
  std::vector<float> pfrechit_z; 
  std::vector<float> pfrechit_x; 
  std::vector<float> pfrechit_y; 
  
  std::vector<float> ecalrechit_eta;
  std::vector<float> ecalrechit_phi;
  std::vector<float> ecalrechit_energy;

  std::vector<float> pcalohit_eta;
  std::vector<float> pcalohit_phi;
  std::vector<float> pcalohit_energy;

  TTree* pfrechit_tree_;
  TTree* ecalrechit_tree_;
  TTree* simhit_tree_;
};

void RecHitDumper::clearVariables() {
  pfrechit_eta.clear();
  pfrechit_phi.clear();
  pfrechit_ieta.clear();
  pfrechit_iphi.clear();
  pfrechit_energy.clear();
  pfrechit_time.clear();
  pfrechit_x.clear();
  pfrechit_y.clear();
  pfrechit_z.clear();

  ecalrechit_eta.clear();
  ecalrechit_phi.clear();
  ecalrechit_energy.clear();

  pcalohit_eta.clear();
  pcalohit_phi.clear();
  pcalohit_energy.clear();
}

RecHitDumper::RecHitDumper(const edm::ParameterSet& ps)
  : pfrechit_token_(consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("pfrechits"))),
    ecalrechit_token_(consumes<EcalRecHitCollection>(ps.getParameter<edm::InputTag>("ecalrechits"))),
    simhit_token_(consumes<edm::PCaloHitContainer>(ps.getParameter<edm::InputTag>("simhits"))),
    geometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()) {
};

RecHitDumper::~RecHitDumper() { clearVariables(); };

void RecHitDumper::beginJob() {
  edm::Service<TFileService> fs;
  pfrechit_tree_ = fs->make<TTree>("pfrechits", "PFRecHits");

  pfrechit_tree_->Branch("pfrechitEta", &pfrechit_eta);
  pfrechit_tree_->Branch("pfrechitPhi", &pfrechit_phi);
  pfrechit_tree_->Branch("pfrechitEnergy", &pfrechit_energy);
  pfrechit_tree_->Branch("pfrechitTime", &pfrechit_time);
  pfrechit_tree_->Branch("pfrechitIEta", &pfrechit_ieta);
  pfrechit_tree_->Branch("pfrechitIPhi", &pfrechit_iphi);
  pfrechit_tree_->Branch("pfrechitX", &pfrechit_x);
  pfrechit_tree_->Branch("pfrechitY", &pfrechit_y);
  pfrechit_tree_->Branch("pfrechitZ", &pfrechit_z);

  ecalrechit_tree_ = fs->make<TTree>("ecalrechits", "ecalrechits");
  ecalrechit_tree_->Branch("ecalrechitEta", &ecalrechit_eta);
  ecalrechit_tree_->Branch("ecalrechitPhi", &ecalrechit_phi);
  ecalrechit_tree_->Branch("ecalrechitEnergy", &ecalrechit_energy);

  simhit_tree_ = fs->make<TTree>("simhits", "simhits");
  simhit_tree_->Branch("simhitEta", &pcalohit_eta);
  simhit_tree_->Branch("simhitPhi", &pcalohit_phi);
  simhit_tree_->Branch("simhitEnergy", &pcalohit_energy);
  
}

void RecHitDumper::beginRun(edm::Run const&, edm::EventSetup const& es) {
  geom = &(es.getData(geometry_token_));
}

void RecHitDumper::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  clearVariables();

  edm::Handle<reco::PFRecHitCollection> pfrechits_h;
  event.getByToken(pfrechit_token_, pfrechits_h);
  const auto& pfrechits = *pfrechits_h;
  
  npfrechits_ = pfrechits.size();
  pfrechit_eta.reserve(npfrechits_);
  pfrechit_phi.reserve(npfrechits_);
  pfrechit_ieta.reserve(npfrechits_);
  pfrechit_iphi.reserve(npfrechits_);
  pfrechit_energy.reserve(npfrechits_);
  pfrechit_time.reserve(npfrechits_);

  std::unordered_map<DetId, float> hitmap;

  for (const auto& hit : pfrechits) {
    uint32_t rawId = hit.detId();
    DetId id(rawId);
    GlobalPoint position = geom->getPosition(id);
    pfrechit_eta.push_back(position.eta());
    pfrechit_phi.push_back(position.phi());
    pfrechit_x.push_back(position.x());
    pfrechit_y.push_back(position.y());
    pfrechit_z.push_back(position.z());

    float energy = hit.energy();
    pfrechit_energy.push_back(energy);
    pfrechit_time.push_back(hit.time());
    hitmap.insert(std::make_pair(id, energy));

    if (id.subdetId() == EcalBarrel) {
      EBDetId eb_id(id);
      pfrechit_ieta.push_back(eb_id.ieta());
      pfrechit_iphi.push_back(eb_id.iphi());
    }
  }
    
  edm::Handle<EcalRecHitCollection> ecalrechits_h;
  event.getByToken(ecalrechit_token_, ecalrechits_h);
  const auto& ecalrechits = *ecalrechits_h;
  for (const auto& hit : ecalrechits) {
    GlobalPoint position = geom->getPosition(hit.detid());
    ecalrechit_eta.push_back(position.eta());
    ecalrechit_phi.push_back(position.phi());
    ecalrechit_energy.push_back(hit.energy());
  }

  edm::Handle<edm::PCaloHitContainer> simhits_h;
  event.getByToken(simhit_token_, simhits_h);
  const auto& simhits = *simhits_h;

  std::unordered_map<uint32_t, double> simhitsMap;
  for (const auto& hit : simhits) {
    std::unordered_map<uint32_t, double>::iterator it = simhitsMap.find(hit.id());
    if (it != simhitsMap.end()) {
      auto energy = it->second + hit.energy();
      it->second = energy;
    } else {
      simhitsMap.emplace(hit.id(), hit.energy());
    }
  }
  for (const auto& hit : simhitsMap) {
    GlobalPoint position = geom->getPosition(DetId(hit.first));
    pcalohit_eta.push_back(position.eta());
    pcalohit_phi.push_back(position.phi());
    pcalohit_energy.push_back(static_cast<float>(hit.second));
  }

  pfrechit_tree_->Fill();
  ecalrechit_tree_->Fill();
  simhit_tree_->Fill();
}

void RecHitDumper::endJob() {}

void RecHitDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pfrechits", edm::InputTag("particleFlowRecHitECAL", ""));
  desc.add<edm::InputTag>("ecalrechits", edm::InputTag("ecalRecHit", "EcalRecHitsEB"));
  desc.add<edm::InputTag>("simhits", edm::InputTag("g4SimHits", "EcalHitsEB"));
  descriptions.add("recHitDumper", desc);
}

DEFINE_FWK_MODULE(RecHitDumper);
  
