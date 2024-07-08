#ifndef __RecoParticleFlow_PFClusterProducer_BarrelLayerClusterProducer_H__
#define __RecoParticleFlow_PFClusterProducer_BarrelLayerClusterProducer_H__

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "RecoParticleFlow/PFClusterProducer/interface/RecHitTopologicalCleanerBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/SeedFinderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterBuilderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFCPositionCalculatorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/CaloRecHitResolutionProvider.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalLayerClusterAlgoFactory.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "RecoParticleFlow/PFClusterProducer/plugins/BarrelCLUEAlgo.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/EBTilesConstants.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HBTilesConstants.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HOTilesConstants.h"

using Density = hgcal_clustering::Density;

class BarrelLayerClusterProducer : public edm::stream::EDProducer<> {
public:
  BarrelLayerClusterProducer(const edm::ParameterSet&);
  ~BarrelLayerClusterProducer() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  reco::CaloCluster::AlgoId algoId;
  std::string timeClname;
  unsigned int nHitsTime;
  edm::EDGetTokenT<reco::PFRecHitCollection> ebhits_token_;
  edm::EDGetTokenT<reco::PFRecHitCollection> hbhits_token_;
  edm::EDGetTokenT<reco::PFRecHitCollection> hohits_token_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

  std::unique_ptr<HGCalClusteringAlgoBase> ebalgo, hbalgo, hoalgo;
  hgcal::RecHitTools rhtools_;
  std::unique_ptr<CaloRecHitResolutionProvider> timeResolutionCalc_;
};

DEFINE_FWK_MODULE(BarrelLayerClusterProducer);

BarrelLayerClusterProducer::BarrelLayerClusterProducer(const edm::ParameterSet& ps)
    : algoId(reco::CaloCluster::undefined),
      timeClname(ps.getParameter<std::string>("timeClname")),
      nHitsTime(ps.getParameter<unsigned int>("nHitsTime")),
      ebhits_token_{consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("EBInput"))},
      hbhits_token_{consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("HBInput"))},
      hohits_token_{consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("HOInput"))},
      caloGeomToken_{consumesCollector().esConsumes<CaloGeometry, CaloGeometryRecord>()} {
  auto ebpluginPSet = ps.getParameter<edm::ParameterSet>("ebplugin");
  ebalgo = HGCalLayerClusterAlgoFactory::get()->create(ebpluginPSet.getParameter<std::string>("type"), ebpluginPSet);
  ebalgo->setAlgoId(algoId);

  auto hbpluginPSet = ps.getParameter<edm::ParameterSet>("hbplugin");
  hbalgo = HGCalLayerClusterAlgoFactory::get()->create(hbpluginPSet.getParameter<std::string>("type"), hbpluginPSet);
  hbalgo->setAlgoId(algoId);

  auto hopluginPSet = ps.getParameter<edm::ParameterSet>("hoplugin");
  hoalgo = HGCalLayerClusterAlgoFactory::get()->create(hopluginPSet.getParameter<std::string>("type"), hopluginPSet);
  hoalgo->setAlgoId(algoId);

  ebalgo->setThresholds(consumesCollector().esConsumes<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd>(),
                        consumesCollector().esConsumes<HcalPFCuts, HcalPFCutsRcd>());
  hbalgo->setThresholds(consumesCollector().esConsumes<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd>(),
                        consumesCollector().esConsumes<HcalPFCuts, HcalPFCutsRcd>());
  hoalgo->setThresholds(consumesCollector().esConsumes<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd>(),
                        consumesCollector().esConsumes<HcalPFCuts, HcalPFCutsRcd>());

  timeResolutionCalc_ = std::make_unique<CaloRecHitResolutionProvider>(ps.getParameterSet("timeResolutionCalc"));
  produces<std::vector<float>>("InitialLayerClustersMaskECAL");
  produces<std::vector<float>>("InitialLayerClustersMaskHCAL");
  produces<std::vector<reco::BasicCluster>>("ecalLayerClusters");
  produces<std::vector<reco::BasicCluster>>("hcalLayerClusters");
  //density
  //produces<Density>();
  //time for layer clusters
  produces<edm::ValueMap<std::pair<float, float>>>(timeClname+"_ecal");
  produces<edm::ValueMap<std::pair<float, float>>>(timeClname+"_hcal");
}

void BarrelLayerClusterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  edm::ParameterSetDescription ebpluginDesc;
  ebpluginDesc.addNode(edm::PluginDescription<HGCalLayerClusterAlgoFactory>("type", "EBCLUE", true));

  edm::ParameterSetDescription hbpluginDesc;
  hbpluginDesc.addNode(edm::PluginDescription<HGCalLayerClusterAlgoFactory>("type", "HBCLUE", true));

  edm::ParameterSetDescription hopluginDesc;
  hopluginDesc.addNode(edm::PluginDescription<HGCalLayerClusterAlgoFactory>("type", "HOCLUE", true));

  desc.add<edm::ParameterSetDescription>("ebplugin", ebpluginDesc);
  desc.add<edm::ParameterSetDescription>("hbplugin", hbpluginDesc);
  desc.add<edm::ParameterSetDescription>("hoplugin", hopluginDesc);

  desc.add<edm::InputTag>("EBInput", edm::InputTag("particleFlowRecHitECAL", ""));
  desc.add<edm::InputTag>("HBInput", edm::InputTag("particleFlowRecHitHBHE", ""));
  desc.add<edm::InputTag>("HOInput", edm::InputTag("particleFlowRecHitHO", ""));

  edm::ParameterSetDescription timeResolutionCalcDesc;
  timeResolutionCalcDesc.addNode(edm::ParameterDescription<double>("noiseTerm", 1.10889, true) and
                                 edm::ParameterDescription<double>("constantTerm", 0.428192, true) and
                                 edm::ParameterDescription<double>("corrTermLowE", 0.0510871, true) and
                                 edm::ParameterDescription<double>("threshLowE", 0.5, true) and
                                 edm::ParameterDescription<double>("constantTermLowE", 0.0, true) and
                                 edm::ParameterDescription<double>("noiseTermLowE", 1.31883, true) and
                                 edm::ParameterDescription<double>("threshHighE", 5.0, true));
  desc.add<edm::ParameterSetDescription>("timeResolutionCalc", timeResolutionCalcDesc);

  desc.add<std::string>("timeClname", "timeLayerCluster");
  desc.add<unsigned int>("nHitsTime", 3);
  descriptions.add("barrelLayerClusters", desc);
}

void BarrelLayerClusterProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  edm::Handle<reco::PFRecHitCollection> ebhits, hbhits, hohits;
  edm::ESHandle<CaloGeometry> geom = es.getHandle(caloGeomToken_);
  rhtools_.setGeometry(*geom);

  /*std::unique_ptr<std::vector<reco::BasicCluster>> ebclusters(new std::vector<reco::BasicCluster>),
						   hbclusters(new std::vector<reco::BasicCluster>),
						   hoclusters(new std::vector<reco::BasicCluster>);
  */
  //std::vector<reco::BasicCluster> ebclusters, hbclusters, hoclusters;

  //auto density = std::make_unique<Density>();
  ebalgo->getEventSetup(es, rhtools_);

  //make a map detid-rechit
  // NB for the moment just host EE and FH hits
  // timing in digi for BH not implemented for now
  std::unordered_map<uint32_t, const reco::PFRecHit*> hitmap;

  evt.getByToken(ebhits_token_, ebhits);
  ebalgo->populate(*ebhits);
  for (auto& hit : *ebhits) {
    hitmap[hit.detId()] = &(hit);
  }
  ebalgo->makeClusters();
  //*ebclusters = ebalgo->getClusters(false);
  std::vector<reco::BasicCluster> ebclusters = ebalgo->getClusters(false);

  hbalgo->getEventSetup(es, rhtools_);

  evt.getByToken(hbhits_token_, hbhits);
  hbalgo->populate(*hbhits);
  for (auto& hit : *hbhits) {
    hitmap[hit.detId()] = &(hit);
  }
  hbalgo->makeClusters();
  //*hbclusters = hbalgo->getClusters(false);
  std::vector<reco::BasicCluster> hbclusters = hbalgo->getClusters(false);
  hoalgo->getEventSetup(es, rhtools_);

  evt.getByToken(hohits_token_, hohits);
  hoalgo->populate(*hohits);
  for (auto& hit : *hohits) {
    hitmap[hit.detId()] = &(hit);
  }
  hoalgo->makeClusters();
  //*hoclusters = hoalgo->getClusters(false);
  std::vector<reco::BasicCluster> hoclusters = hoalgo->getClusters(false);
  
  std::unique_ptr<std::vector<reco::BasicCluster>> clusters_ecal(new std::vector<reco::BasicCluster>);
  (*clusters_ecal).reserve(ebclusters.size());
  (*clusters_ecal).insert((*clusters_ecal).end(), ebclusters.begin(), ebclusters.end());

  std::unique_ptr<std::vector<reco::BasicCluster>> clusters_hcal(new std::vector<reco::BasicCluster>);
  (*clusters_hcal).reserve(hbclusters.size());
  (*clusters_hcal).insert((*clusters_hcal).end(), hbclusters.begin(), hbclusters.end());
  //(*clusters).insert((*clusters).end(), hoclusters.begin(), hoclusters.end());
  
  auto clusterHandleEcal = evt.put(std::move(clusters_ecal), "ecalLayerClusters");
  auto clusterHandleHcal = evt.put(std::move(clusters_hcal), "hcalLayerClusters");
  //Keep the density
  /**density = algo->getDensity();
  evt.put(std::move(density));*/

  edm::PtrVector<reco::BasicCluster> clusterEcalPtrs;  //, clusterPtrsSharing;
  edm::PtrVector<reco::BasicCluster> clusterHcalPtrs;  //, clusterPtrsSharing;
  
  std::vector<std::pair<float, float>> times_ecal, times_hcal;
  times_ecal.reserve(clusterHandleEcal->size());
  times_hcal.reserve(clusterHandleHcal->size());

  for (unsigned i = 0; i < clusterHandleEcal->size(); ++i) {
    edm::Ptr<reco::BasicCluster> ptr(clusterHandleEcal, i);
    clusterEcalPtrs.push_back(ptr);

    std::pair<float, float> timeCl(-99., -1.);

    const reco::CaloCluster& sCl = (*clusterHandleEcal)[i];
    if (sCl.size() >= nHitsTime) {
      std::vector<float> timeClhits;
      std::vector<float> timeErrorClhits;

      for (auto const& hit : sCl.hitsAndFractions()) {
        auto finder = hitmap.find(hit.first);
        if (finder == hitmap.end())
          continue;

        //time is computed wrt  0-25ns + offset and set to -1 if no time
        const reco::PFRecHit* rechit = finder->second;
        float rhTimeE = timeResolutionCalc_->timeResolution2(rechit->energy());
        //check on timeError to exclude scintillator
        if (rhTimeE < 0.)
          continue;
        timeClhits.push_back(rechit->time());
        timeErrorClhits.push_back(1. / rhTimeE);
      }
      hgcalsimclustertime::ComputeClusterTime timeEstimator;
      timeCl = timeEstimator.fixSizeHighestDensity(timeClhits, timeErrorClhits, nHitsTime);
    }
    times_ecal.push_back(timeCl);
  }

  for (unsigned i = 0; i < clusterHandleHcal->size(); ++i) {
    edm::Ptr<reco::BasicCluster> ptr(clusterHandleHcal, i);
    clusterHcalPtrs.push_back(ptr);

    std::pair<float, float> timeCl(-99., -1.);

    const reco::CaloCluster& sCl = (*clusterHandleHcal)[i];
    if (sCl.size() >= nHitsTime) {
      std::vector<float> timeClhits;
      std::vector<float> timeErrorClhits;

      for (auto const& hit : sCl.hitsAndFractions()) {
        auto finder = hitmap.find(hit.first);
        if (finder == hitmap.end())
          continue;

        //time is computed wrt  0-25ns + offset and set to -1 if no time
        const reco::PFRecHit* rechit = finder->second;
        float rhTimeE = timeResolutionCalc_->timeResolution2(rechit->energy());
        //check on timeError to exclude scintillator
        if (rhTimeE < 0.)
          continue;
        timeClhits.push_back(rechit->time());
        timeErrorClhits.push_back(1. / rhTimeE);
      }
      hgcalsimclustertime::ComputeClusterTime timeEstimator;
      timeCl = timeEstimator.fixSizeHighestDensity(timeClhits, timeErrorClhits, nHitsTime);
    }
    times_hcal.push_back(timeCl);
  }

  std::unique_ptr<std::vector<float>> layerClustersMaskEcal(new std::vector<float>);
  layerClustersMaskEcal->resize(clusterHandleEcal->size(), 1.0);
  evt.put(std::move(layerClustersMaskEcal), "InitialLayerClustersMaskECAL");

  std::unique_ptr<std::vector<float>> layerClustersMaskHcal(new std::vector<float>);
  layerClustersMaskHcal->resize(clusterHandleHcal->size(), 1.0);
  evt.put(std::move(layerClustersMaskHcal), "InitialLayerClustersMaskHCAL");
    
  auto timeCl_ecal = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
  edm::ValueMap<std::pair<float, float>>::Filler filler_ecal(*timeCl_ecal);
  filler_ecal.insert(clusterHandleEcal, times_ecal.begin(), times_ecal.end());
  filler_ecal.fill();
  evt.put(std::move(timeCl_ecal), timeClname+"_ecal");

   auto timeCl_hcal = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
   edm::ValueMap<std::pair<float, float>>::Filler filler_hcal(*timeCl_hcal);
   filler_hcal.insert(clusterHandleHcal, times_hcal.begin(), times_hcal.end());
   filler_hcal.fill();
   evt.put(std::move(timeCl_hcal), timeClname+"_hcal");
   
  /*if (doSharing) {
    for (unsigned i = 0; i < clusterHandleSharing->size(); ++i) {
      edm::Ptr<reco::BasicCluster> ptr(clusterHandleSharing, i);
      clusterPtrsSharing.push_back(ptr);
    }
  }*/
  ebalgo->reset();
  hbalgo->reset();
  hoalgo->reset();
}

#endif  //__RecoLocalCalo_HGCRecProducers_BarrelLayerClusterProducer_H__
