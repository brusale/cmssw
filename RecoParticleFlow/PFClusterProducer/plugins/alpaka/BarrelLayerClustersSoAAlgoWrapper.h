#ifndef RecoParticleFlow_PFClusterProducer_plugins_alpaka_BarrelLayerClustersSoAAlgoWrapper_h
#define RecoParticleFlow_PFClusterProducer_plugins_alpaka_BarrelLayerClustersSoAAlgoWrapper_h

#include <Eigen/Core>
#include <alpaka/alpaka.hpp>

#include "DataFormats/ParticleFlowReco/interface/PFRecHitHostCollection.h"
#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitDeviceCollection.h"
#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitExtraDeviceCollection.h"
#include "DataFormats/ParticleFlowReco/interface/alpaka/BarrelSoAClustersDeviceCollection.h"
#include "RecoParticleFlow/PFClusterProducer/interface/alpaka/BarrelSoAClustersExtraDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class BarrelLayerClustersSoAAlgoWrapper {
  public:
    void run(Queue& queue,
             const unsigned int numer_of_clusters,
             float thresholdW0,
             float positionDeltaRho2,
             const reco::PFRecHitDeviceCollection::ConstView input_rechits_soa,
             const reco::PFRecHitExtraDeviceCollection::ConstView input_clusters_soa,
             BarrelSoAClustersDeviceCollection::View outputs,
             BarrelSoAClustersExtraDeviceCollection::View outputs_service) const;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoParticleFlow_PFClusterProducer_plugins_alpaka_BarrelLayerClustersSoAAlgoWrapper_h
