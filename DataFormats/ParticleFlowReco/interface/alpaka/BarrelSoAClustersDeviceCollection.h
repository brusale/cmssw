#ifndef DataFormats_PortableTestObjects_interface_alpaka_BarrelSoAClustersDeviceCollection_h
#define DataFormats_PortableTestObjects_interface_alpaka_BarrelSoAClustersDeviceCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/ParticleFlowReco/interface/BarrelSoAClusters.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using BarrelSoAClustersDeviceCollection = PortableCollection<BarrelSoAClusters>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // DataFormats_PortableTestObjects_interface_alpaka_BarrelSoAClustersDeviceCollection_h
