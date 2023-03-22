// user include files
#include <unordered_map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

class BarrelRecHitMapProducer : public edm::global::EDProducer<> {
public:
  BarrelRecHitMapProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  const edm::EDGetTokenT<reco::PFRecHitCollection> hits_eb_token_;
  const edm::EDGetTokenT<reco::PFRecHitCollection> hits_hb_token_;
  const edm::EDGetTokenT<reco::PFRecHitCollection> hits_ho_token_;
};

DEFINE_FWK_MODULE(BarrelRecHitMapProducer);

using DetIdRecHitMap = std::unordered_map<DetId, const reco::PFRecHit*>;

BarrelRecHitMapProducer::BarrelRecHitMapProducer(const edm::ParameterSet& ps)
    : hits_eb_token_(consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("EBInput"))),
      hits_hb_token_(consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("HBInput"))),
      hits_ho_token_(consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("HOInput"))) {
  produces<DetIdRecHitMap>();
}

void BarrelRecHitMapProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("EBInput", {"particleFlowRecHitECAL", ""});
  desc.add<edm::InputTag>("HBInput", {"particleFlowRecHitHBHE", ""});
  desc.add<edm::InputTag>("HOInput", {"particleFlowRecHitHO", ""});
  descriptions.add("barrelRecHitMapProducer", desc);
}

void BarrelRecHitMapProducer::produce(edm::StreamID, edm::Event& evt, const edm::EventSetup& es) const {
  auto hitMap = std::make_unique<DetIdRecHitMap>();
  const auto& eb_hits = evt.get(hits_eb_token_);
  const auto& hb_hits = evt.get(hits_hb_token_);
  const auto& ho_hits = evt.get(hits_ho_token_);

  for (const auto& hit : eb_hits) {
    hitMap->emplace(hit.detId(), &hit);
  }

  for (const auto& hit : hb_hits) {
    hitMap->emplace(hit.detId(), &hit);
  }

  for (const auto& hit : ho_hits) {
    hitMap->emplace(hit.detId(), &hit);
  }
  evt.put(std::move(hitMap));
}
