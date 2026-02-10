#include "DataFormats/ParticleFlowReco/interface/BarrelSoAClusters.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitHostCollection.h"
#include "DataFormats/ParticleFlowReco/interface/alpaka/BarrelSoAClustersDeviceCollection.h"
#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitExtraDeviceCollection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoParticleFlow/PFClusterProducer/interface/BarrelSoAClustersExtra.h"
#include "RecoParticleFlow/PFClusterProducer/interface/alpaka/BarrelSoAClustersExtraDeviceCollection.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/BarrelTilesConstants.h"

#include "BarrelLayerClustersSoAAlgoWrapper.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class BarrelSoALayerClustersProducer : public stream::SynchronizingEDProducer<> {
    public:
      BarrelSoALayerClustersProducer(edm::ParameterSet const& config)
        : SynchronizingEDProducer(config),
          getTokenDeviceRecHits_{consumes(config.getParameter<edm::InputTag>("pfRecHitsSoA"))},
          getTokenDeviceClusters_{consumes(config.getParameter<edm::InputTag>("pfRecHitsLayerClustersSoA"))},
          deviceTokenSoAClusters_{produces()} {}

      ~BarrelSoALayerClustersProducer() override = default;

      void acquire(device::Event const& iEvent, device::EventSetup const& iSetup) override {
        auto const& deviceInputClusters = iEvent.get(getTokenDeviceClusters_);
        auto const inputClusters_v = deviceInputClusters.view();

        auto device_numclusters = cms::alpakatools::make_device_view<const unsigned int>(
          alpaka::getDev(iEvent.queue()), inputClusters_v.numberOfClustersScalar());
        auto host_numclusters = cms::alpakatools::make_host_view<unsigned int>(num_clusters_);
        alpaka::memcpy(iEvent.queue(), host_numclusters, device_numclusters);
      }

      void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
        printf("%s %i\n", __FILE__, __LINE__);
        auto const& deviceInputRecHits = iEvent.get(getTokenDeviceRecHits_);
        printf("%s %i\n", __FILE__, __LINE__);
        auto const& inputRechits_v = deviceInputRecHits.view();
        printf("%s %i\n", __FILE__, __LINE__);

        auto const& deviceInputClusters = iEvent.get(getTokenDeviceClusters_);
        printf("%s %i\n", __FILE__, __LINE__);
        auto const& inputClusters_v = deviceInputClusters.view();
        printf("%s %i\n", __FILE__, __LINE__);

        BarrelSoAClustersDeviceCollection output(num_clusters_, iEvent.queue());
        printf("%s %i\n", __FILE__, __LINE__);
        auto output_v = output.view();
        printf("%s %i\n", __FILE__, __LINE__);
        BarrelSoAClustersExtraDeviceCollection outputWorkspace(num_clusters_, iEvent.queue());
        printf("%s %i\n", __FILE__, __LINE__);
        auto output_workspace_v = outputWorkspace.view();
        printf("%s %i\n", __FILE__, __LINE__);

        algo_.run(iEvent.queue(),
                  num_clusters_,
                  0, // thresholdW0 in HGCAL
                  999999, // positionDeltaRho2 in HGCAL
                  inputRechits_v,
                  inputClusters_v,
                  output_v,
                  output_workspace_v);
        printf("%s %i\n", __FILE__, __LINE__);
        iEvent.emplace(deviceTokenSoAClusters_, std::move(output));
        printf("%s %i\n", __FILE__, __LINE__);
      }

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
        edm::ParameterSetDescription desc;
        desc.add<edm::InputTag>("pfRecHitsLayerClustersSoA", edm::InputTag("barrelSoARecHitsLayerClustersProducer"));
        desc.add<edm::InputTag>("pfRecHitsSoA", edm::InputTag("pfRecHitsSoA"));
        descriptions.addWithDefaultLabel(desc);
      }
    
    private:
      device::EDGetToken<reco::PFRecHitDeviceCollection> const getTokenDeviceRecHits_;
      device::EDGetToken<reco::PFRecHitExtraDeviceCollection> const getTokenDeviceClusters_;
      device::EDPutToken<BarrelSoAClustersDeviceCollection> const deviceTokenSoAClusters_;
      BarrelLayerClustersSoAAlgoWrapper algo_;
      unsigned int num_clusters_;
  };
} // namespace ALPAKA_ACCELERATOR_NAMESPACE


#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(BarrelSoALayerClustersProducer);
