#ifndef DataFormats_PortableTestObjects_interface_alpaka_PFRecHitExtraDeviceCollection_h
#define DataFormats_PortableTestObjects_interface_alpaka_PFRecHitExtraDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitExtra.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitExtraHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::reco {

  using ::reco::PFRecHitExtraHostCollection;

  // SoA with CLUE algo parameters (delta, rho, isSeed, ...)
  using PFRecHitExtraDeviceCollection = PortableCollection<::reco::PFRecHitExtra>;
} // namespace ALPAKA_ACCELERATOR_NAMESPACE

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(reco::PFRecHitExtraDeviceCollection, reco::PFRecHitExtraHostCollection);

#endif
