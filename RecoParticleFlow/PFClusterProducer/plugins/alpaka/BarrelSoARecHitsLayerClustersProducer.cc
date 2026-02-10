#include "DataFormats/ParticleFlowReco/interface/PFRecHitHostCollection.h"
#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitDeviceCollection.h"
#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitExtraDeviceCollection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/BarrelTilesConstants.h"

#include "BarrelLayerClustersAlgoWrapper.h"

// Processes the input RecHit SoA collection and generates an output SoA
// containing all the necessary information to build the clusters.
// Specifically, this producer does not create the clusters in any format.
// Instead, it fills a SoA (HGCalSoARecHitsExtra) with the same size as the input
// RecHit SoA. This output SoA includes all the data needed to assemble the
// clusters and assigns a clusterId to each cell that belongs to a cluster.
// Consequently, this producer must be used by another downstream producer to
// either build traditional clusters or to create a SoA representing the
// clusters, complete with all required information (e.g., energy, position).
namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class BarrelSoARecHitsLayerClustersProducer : public stream::EDProducer<> {
  public:
    BarrelSoARecHitsLayerClustersProducer(edm::ParameterSet const& config)
        : EDProducer(config),
          getTokenDevice_{consumes(config.getParameter<edm::InputTag>("pfRecHitsSoA"))},
          deviceToken_{produces()},
          deltac_((float)config.getParameter<double>("deltac")),
          kappa_((float)config.getParameter<double>("kappa")),
          outlierDeltaFactor_((float)config.getParameter<double>("outlierDeltaFactor")) {}

    ~BarrelSoARecHitsLayerClustersProducer() override = default;

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      auto const& deviceInput = iEvent.get(getTokenDevice_);
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      //std::cout << "Size of device collection: " << deviceInput->metadata().size() << std::endl;
      auto const input_v = deviceInput.view();
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      // Allocate output SoA
      reco::PFRecHitExtraDeviceCollection output(deviceInput->metadata().size(), iEvent.queue());
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      auto output_v = output.view();
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      std::cout << "deviceInput->metadata().size(): " << deviceInput->metadata().size() << std::endl;
      std::cout << "deltac_: " << deltac_ << std::endl;
      std::cout << "kappa_: " << kappa_ << std::endl;
      std::cout << "outlierDeltaFactor_: " << outlierDeltaFactor_ << std::endl;
      algo_.run(
          iEvent.queue(), deviceInput->size(), deltac_, kappa_, outlierDeltaFactor_, input_v, output_v);
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      iEvent.emplace(deviceToken_, std::move(output));
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("pfRecHitsSoA", edm::InputTag("TO BE DEFINED"));
      desc.add<double>("deltac", 1.3);
      desc.add<double>("kappa", 9.);
      desc.add<double>("outlierDeltaFactor", 2.);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    // use device::EDGetToken<T> to read from device memory space
    device::EDGetToken<reco::PFRecHitDeviceCollection> const getTokenDevice_;
    device::EDPutToken<reco::PFRecHitExtraDeviceCollection> const deviceToken_;
    BarrelLayerClustersAlgoWrapper algo_;
    const float deltac_;
    const float kappa_;
    const float outlierDeltaFactor_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(BarrelSoARecHitsLayerClustersProducer);
