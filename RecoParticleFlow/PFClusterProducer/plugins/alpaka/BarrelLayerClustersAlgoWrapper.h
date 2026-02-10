#ifndef RecoParticleFlow_PFClusterProducer_plugins_alpaka_BarrelLayerClustersAlgoWrapper_h
#define RecoParticleFlow_PFClusterProducer_plugins_alpaka_BarrelLayerClustersAlgoWrapper_h

#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitDeviceCollection.h"
#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitExtraDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class BarrelLayerClustersAlgoWrapper {
  public:
    void run(Queue& queue,
             const unsigned int size,
             const float dc,
             const float kappa,
             const float outlierDeltaFactor,
             const reco::PFRecHitDeviceCollection::ConstView inputs,
             reco::PFRecHitExtraDeviceCollection::View outputs) const;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoParticleFlow_PFClusterProducer_plugins_alpaka_HGCalLayerClustersAlgoWrapper_h
