// Original author: Leonardo Cristella

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociator.h"
#include "SimBarrelLCToSCAssociatorByEnergyScoreImpl.h"

class SimBarrelLCToSCAssociatorByEnergyScoreProducer : public edm::global::EDProducer<> {
public:
  explicit SimBarrelLCToSCAssociatorByEnergyScoreProducer(const edm::ParameterSet &);
  ~SimBarrelLCToSCAssociatorByEnergyScoreProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  edm::EDGetTokenT<std::unordered_map<DetId, float>> hitMap_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometry_;
  const bool hardScatterOnly_;
  std::shared_ptr<hgcal::RecHitTools> rhtools_;
};

SimBarrelLCToSCAssociatorByEnergyScoreProducer::SimBarrelLCToSCAssociatorByEnergyScoreProducer(const edm::ParameterSet &ps)
    : hitMap_(consumes<std::unordered_map<DetId, float>>(ps.getParameter<edm::InputTag>("hitMapTag"))),
      caloGeometry_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      hardScatterOnly_(ps.getParameter<bool>("hardScatterOnly")) {
  rhtools_.reset(new hgcal::RecHitTools());

  // Register the product
  produces<hgcal::LayerClusterToSimClusterAssociator>();
}

SimBarrelLCToSCAssociatorByEnergyScoreProducer::~SimBarrelLCToSCAssociatorByEnergyScoreProducer() {}

void SimBarrelLCToSCAssociatorByEnergyScoreProducer::produce(edm::StreamID,
                                                    edm::Event &iEvent,
                                                    const edm::EventSetup &es) const {
  edm::ESHandle<CaloGeometry> geom = es.getHandle(caloGeometry_);
  rhtools_->setGeometry(*geom);

  const auto hitMap = &iEvent.get(hitMap_);

  auto impl =
      std::make_unique<SimBarrelLCToSCAssociatorByEnergyScoreImpl>(iEvent.productGetter(), hardScatterOnly_, rhtools_, hitMap);
  auto toPut = std::make_unique<hgcal::LayerClusterToSimClusterAssociator>(std::move(impl));
  iEvent.put(std::move(toPut));
}

void SimBarrelLCToSCAssociatorByEnergyScoreProducer::fillDescriptions(edm::ConfigurationDescriptions &cfg) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("hitMapTag", edm::InputTag("simBarrelRecHitMapProducer"));
  desc.add<bool>("hardScatterOnly", true);

  cfg.add("simBarrelSimClusterAssociatorByEnergyScore", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimBarrelLCToSCAssociatorByEnergyScoreProducer);
