#ifndef DataFormats_CaloRecHit_CaloClustersDeviceCollection_h
#define DataFormats_CaloRecHit_CaloClustersDeviceCollection_h

#include "DataFormats/CaloRecHit/interface/CaloClusters.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using CaloClustersDeviceCollection = PortableCollection<reco::CaloCluster>;
}

#endif
