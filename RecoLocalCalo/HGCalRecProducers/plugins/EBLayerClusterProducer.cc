#ifndef __RecoLocalCalo_HGCalRecProducers_EBLayerClusterProducer_H__
#define __RecoLocalCalo_HGCalRecProducers_EBLayerClusterProducer_H__

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescription.h"
#include "FWCore/ParameterSet/inteface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"


#include "RecoLocalCalo/HGCalRecProducer/interface/HGCalLayerClusterAlgoFactory.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/Common/interface/ValueMap.h"

using Density = hgcal_clustering::Density;

class EBLayerClusterProducer : public edm::stream::EDProducer<> {
  public:
    EBLayerClusterProducer(const edm::ParameterSet&);
    ~EBLayerClusterProducer() override {}
    static void fillDescriptions(edm::ConfigurationDescription& descriptions);
  
    void produce(edm::Event&, const edm::EventSetup&) override;

  private:
    edm::EDGetTokenT<reco::PFRecHitCollection> hits_token_;

    std::unique_ptr<HGCalClusteringAlgoBase> algo_;
    std::string timeClName_; 
    unsigned int nHitsTime_;

DEFINE_FWK_MODULE(EBLayerClusterProducer);

EBLayerClusterProducer::EBLayerClusterProducer(const edm::ParameterSet& ps)
  : timeClName_(ps.getParameter<std::string>("timeClName")),
    nHitsTime_(ps.getParameter<unsigned int>("nHitsTime")) {
  hits_token_ = consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("EBInput"));

  auto pluginPSet = ps.getParameter<edm::ParameterSet>("plugin");
  algo_ = HGCalLayerClusterAlgoFactory::get->create(
    pluginPSet.getParameter<std::string>("type"), pluginPSet, consumesCollector());

  produces<std::vector<float>>("InitialLayerClustersMask");
  produces<std::vector<reco::BasicCluster>>();
  produces<Density>();
  produces<edm::ValueMap<std::pair<float, float>>>(timeClName);

}

void EBLayerClusterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  edm::ParameterSetDescription pluginDesc;
  pluginDesc.addNode(edm::PluginDescription<HGCalLayerClusterAlgoFactory>("type", "CLUE", true));

  desc.add<edm::ParameterSetDescription>("plugin", pluginDesc);
  desc.add<edm::InputTag>("EBInput", edm::InputTag("particleFlowRecHitECAL", "Cleaned"));
  desc.add<std::string>("timeClName", "timeLayerCluster");
  desc.add<unsigned int>("nHitsTime", 3);

  descriptions.add("ebLayerClusters", desc);
}

void EBLayerClusterProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  edm::Handle<reco::PFRecHitCollection> hits;

  std::unique_ptr<std::vector<reco::BasicCluster>> clusters(new std::vector<reco::BasicCluster>);
  auto density = std::make_unique<Density>;

  //algo->getEventSetup();

  std::unordered_map<uint32_t, const PFRecHit*> hitmap;

  evt.getByToken(hits_token, hits);
  
  algo_->populate(hits);
  for (auto hit : hits) {
    hitmap[hit.detid()] = &(hit);
  }
  algo_->makeClusters();
  *clusters = algo_->getClusters(false);

  auto clusterHandle = evt.put(std::move(clusters));

  *density = algo_->getDensity();
  evt.put(std::move(density));

  edm::PtrVector<reco::BasicCluster> clusterPtrs;

  std::vector<std::pair<float, float> times;
  times.reserve(clusterHandle->size());

  for (unsigned i = 0; i < clusterHandle->size(); ++i) {
    edm::Ptr<reco::BasicCluster> ptr(clusterHandle, i);
    clusterPtrs.push_back(ptr);

    std::pair<float, float> timeCl(-99., -1.);

    const reco::CaloCluster& sCl = (*clusterHandle)[i];
    if (sCl.size() >= nHitsTime) {
      std::vector<float> timeClHits;
      std::vector<float> timeErrorClhits;

      for (auto const& hit : sCl.hitsAndFractions()) {
        auto finder = hitmap.find(hit.first);
        if (finder == hitmap.end())
          continue;

        const PFRecHit* rechit = finder->second;
        //float rhTimeE = rechit->timeError(); //understand this, PFRecHit does not have timeError
        if (rhTime < 0) 
          continue;
        timeClHits.push_back(rechit->time());
        //timeErrorClhits.push_back(1./(rhTimeE*rhtimeE));
      }
      hgcalsimclustertime::ComputeClusterTime timeEstimator;
      timeCl = timeEstimator.fixSizeHighestDensity(timeClhits, timerErrorClhits, nHitsTime);
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
  evt.put(std::move(timeCl), timeClname);

  algo->reset();
}

#endif
