// user include files
#include <unordered_map>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

class SimBarrelHitMapProducer : public edm::global::EDProducer<> {
public:
  SimBarrelHitMapProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  const edm::EDGetTokenT<edm::PCaloHitContainer> hits_token_;
};

DEFINE_FWK_MODULE(SimBarrelHitMapProducer);

using DetIdRecHitMap = std::unordered_map<DetId, float>;

SimBarrelHitMapProducer::SimBarrelHitMapProducer(const edm::ParameterSet& ps)
    : hits_token_(consumes<edm::PCaloHitContainer>(ps.getParameter<edm::InputTag>("pCaloHits"))) {
  produces<DetIdRecHitMap>();
}

// maybe take SimHits straight from CP? 
void SimBarrelHitMapProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pCaloHits", {"g4SimHits", "EcalHitsEB"});
  descriptions.add("simBarrelRecHitMapProducer", desc);
}

void SimBarrelHitMapProducer::produce(edm::StreamID, edm::Event& evt, const edm::EventSetup& es) const {
  auto hitMap = std::make_unique<DetIdRecHitMap>();
  const auto& hits = evt.get(hits_token_);

  std::unordered_map<uint32_t, PCaloHit> hit_and_energies;
  for (unsigned int i = 0; i < hits.size(); ++i) {
    PCaloHit hit = hits[i];
    std::unordered_map<uint32_t, PCaloHit>::iterator it = hit_and_energies.find(hit.id());
    if (it != hit_and_energies.end()) {
      double energy = hit.energy() + it->second.energy();
      it->second.setEnergy(energy);
    } else {
      hit_and_energies.emplace(hit.id(), hit);
    }
  }
  
  /*for (auto it = hit_and_energies.begin(); it != hit_and_energies.end(); ++it) {
    if (it->second.energy() > 0.1) {
      std::cout << "New hit found with id: " << it->first << std::endl;
      std::cout << "New hit found with energy: " << it->second.energy() << std::endl;
      DetId id(it->first);
      const PCaloHit hit = it->second;
      DetIdRecHitMap::const_iterator hitmap_it = hitMap->find(id);
      if (hitmap_it == hitMap->end()) {
	hitMap->try_emplace(id, &hit);
	std::cout << "ptr detid: " << id.rawId() << std::endl;
	std::cout << "ptr energy: " << hitMap->at(id)->energy() << std::endl;
      }
    }
  } */

  for (const auto& [key, cH] : hit_and_energies) {
    //if (cH.energy() > 0.1 && (DetId(cH.id()).det() == 3 || DetId(cH.id()).det() == 4)) {
    if ((DetId(cH.id()).det() == 3 || DetId(cH.id()).det() == 4)) {  
      //std::cout << "New hit found with id: " << cH.id() << std::endl;
      //std::cout << "New hit found with energy: " << cH.energy() << std::endl;
      DetId id(cH.id());
      const PCaloHit* hit = &cH;
      hitMap->emplace(id, static_cast<float>(hit->energy()));
      //std::cout << "ptr detid: " << id.rawId() << std::endl;
      //std::cout << "ptr energy: " << hitMap->at(id) << std::endl;
    }
  }

  /*for (auto it = hitMap->begin(); it != hitMap->end(); ++it) {
    std::cout << "Det: " << (it->first).rawId() << std::endl;
    std::cout << "Energy: " << hitMap->at(it->first) << std::endl;
  }*/


  evt.put(std::move(hitMap));
}
