#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include "RecoLocalCalo/HGCalRecProducers/interface/BarrelTilesConstants.h"
#include "BarrelLayerClustersAlgoWrapper.h"

#include "RecoParticleFlow/PFClusterProducer/interface/alpaka/CLUEAlgoAlpaka.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;
  
  void BarrelLayerClustersAlgoWrapper::run(Queue& queue,
                                           const unsigned int size,
                                           const float dc,
                                           const float kappa,
                                           const float outlierDeltaFactor,
                                           const reco::PFRecHitDeviceCollection::ConstView inputs,
                                           reco::PFRecHitExtraDeviceCollection::View outputs) const {
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    CLUEAlgoAlpaka<ALPAKA_ACCELERATOR_NAMESPACE::Acc1D, Queue, EBTilesConstants, 1> algoStandalone(
      queue, dc, kappa, outlierDeltaFactor, true);

    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    algoStandalone.makeClustersCMSSW(size,
                                     inputs.eta().data(),
                                     inputs.phi().data(),
                                     inputs.depth().data(),
                                     inputs.energy().data(),
                                     inputs.sigmaNoise().data(),
                                     inputs.detId().data(),
                                     outputs.rho().data(),
                                     outputs.delta().data(),
                                     outputs.nearestHigher().data(),
                                     outputs.clusterIndex().data(),
                                     outputs.isSeed().data(),
                                     &outputs.numberOfClustersScalar());
  }
} // namespace ALPAKA_ACCELERATOR_NAMESPACE
