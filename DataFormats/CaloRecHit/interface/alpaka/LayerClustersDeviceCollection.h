#ifndef DataFormats_CaloRecHit_LayerClustersDeviceCollection_h
#define DataFormats_CaloRecHit_LayerClustersDeviceCollection_h

#include "DataFormats/CaloRecHit/interface/LayerClusters.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using LayerClustersDeviceCollection = PortableCollection<reco::LayerClusters>;
}

#endif
