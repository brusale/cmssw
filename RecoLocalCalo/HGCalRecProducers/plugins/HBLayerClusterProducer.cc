#ifndef __RecoLocalCalo_HGCalRecProducers_HBLayerClustersProducer_H__
#define __RecoLocalCalo_HGCalRecProducers_HBLayerClusterProducer_H__

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalLayerClusterAlgoFactory.h"
#include "RecoParticleFlow/PFClusterProducer/interface/CaloRecHitResolutionProvider.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"


#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/Common/interface/ValueMap.h"

using Density = hgcal_clustering::Density;

class HBLayerClusterProducer : public edm::stream::EDProducer<> {
  public:
    HBLayerClusterProducer(const edm::ParameterSet&);
    ~HBLayerClusterProducer() override {}
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    void produce(edm::Event&, const edm::EventSetup&) override;

  private:
    edm::EDGetTokenT<reco::PFRecHitCollection> hits_token_;
    edm::EDGetTokenT<reco::PFRecHitCollection> hits_token_;
    std::unique_ptr<HGCalClusteringAlgoBase> algo_;

    std::string timeClName_;
    unsigned int nHitsTime_;
};

DEFINE_FWK_MODULE(HBLayerClusterProducer);

HBLayerClusterProducer::HBLayerClusterProducer(const edm::ParameterSet& ps)
  : timeClName_(ps.getParameter<std::string>("timeClName")),
    nHitsTime_(ps.getParameter<unsigned int>("nHitsTime")) {
  
  hbhits_token_ = consumes<reco::PFRecHitCollection>(ps.getParameter("HBInput"));
  hohits_token_ = consumes<reco::PFRecHitCollection>(ps.getParameter("HOInput"));

  auto pluginPSet = ps.getParameter<edm::ParameterSet>("plugin");
  algo_ = HGCalLayerClusterAlgoFactory::get()->create(
    pluginPSet.getParameter<std::string>("type"), pluginPSet, consumesCollector());

  produces<std::vector<reco::BasicCluster>>();
  //produces<std::vector<float>>("InitialLayerClusterMask");
  //produce<Density>();
  //produces<edm::ValueMap<std::pair<float, float>>>(timeClName_);
}

void HBLayerClusterProducer::fillDescriptions(edm::Configurations& descriptions) {
  edm::ParameterSetDescription desc;
  edm::ParameterSetDescription pluginDesc;
  //edm::ParameterSetDescription timeResDesc;
  pluginDesc.addNode(edm::pluginPSetDescription<HGCalLayerClusterAlgoFactor>("type", "HBCLUE", true));

  desc.add<edm::ParameterSetDescription>("plugin", pluginDesc);
  //desc.add<edm::ParameterSetDescription>("timeDescription", timeResDesc);
  desc.add<edm::InputTag>("HBInput", edm::InputTag("particleFlowRecHitHBHE", "Cleaned"));
  desc.add<std::string>("timeClName", "timeLayerCluster");
  desc.add<unsigned int>("nHitsTime", 3);

  descriptions.add("hbLayerClusters", desc);
}

void HBLayerClusterProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  edm::Handle<reco::PFRecHitCollection> hits;
  std::unique_ptr<std::vector<reco::BasicCluster>> clusters(new std::vector<reco::BasicCluster>);
  //auto density = std::make_unique<Density>();

  algo_->getEventSetup(es);

  std::unordered_map<uint32_t, const reco::PFRecHit*> hitmap;
  evt.getByToken(hits_token_, hits);

  algo_->populate(*hits);
  
  for (auto const& hit : hits) {
    hitmap[hit.detId()] = &(hit);
  }
  algo_->makeClusters();
  *clusters = algo_->getClusters(false);

  auto clusterHandle = evt.put(std::move(clusters));

  *density = algo_->getDensity();
  evt.put(std::move(density));

  edm::PtrVector<reco::BasicCluster> clusterPtrs;

  std::vector<std::pair<float, float>> times;
  times.reserve(clusterHandle->size());

  for (unsigned i = 0; i < clusterHandle->size(); ++i) {
    edm::Ptr<reco::BasicCluster> ptr(clusterHandle, i);
    clusterPtrs.push_back(ptr);

    std::pair<float, float> timeCl(-99., -1.);

    const reco::CaloCluster& sCl = (*clusterHandle)[i];
    if (sCl.size() >= nHitsTime_) {
      std::vector<float> timeClhits;
      std::vector<float> timeErrorClhits;

      for (auto const& hit : sCl.hitsAndFractions()) {
        auto finder = hitmap.find(hit.first);
        if (finder == hitmap.end())
          continue;

        const reco::PFRecHit* rechit = finder->second;
        //float rhTimeE = rechit->timeError(); //understand this, PFRecHit does not have timeError
        //float rhTimeE = sqrt(timeResolutionCalc_->timeResolution2(rechit->energy()));
        //float rhTimeE = rechit->time()/10.; //just for testing purposes
        if (rhTimeE < 0) 
          continue;
        timeClhits.push_back(rechit->time());
        timeErrorClhits.push_back(1./(rhTimeE*rhTimeE));
      }
      hgcalsimclustertime::ComputeClusterTime timeEstimator;
      timeCl = timeEstimator.fixSizeHighestDensity(timeClhits, timeErrorClhits, nHitsTime_);
    }
    times.push_back(timeCl);
  }
  std::unique_ptr<std::vector<float>> layerClustersMask(new std::vector<float>);
  layerClustersMask->resize(clusterHandle->size(), 1.0);
  evt.put(std::move(layerClustersMask), "InitialLayerClustersMaks");

  auto timeCl = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
  edm::ValueMap<std::pair<float, float>>::Filler filler(*timeCl);
  filler.insert(clusterHandle, times.begin(), times.end());
  filler.fill();
  evt.put(std::move(timeCl), timeClName_);


}








#endif
