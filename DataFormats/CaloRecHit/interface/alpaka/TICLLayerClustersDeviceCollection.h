#ifndef DataFormats_CaloRecHit_TICLLayerClustersDeviceCollection_h
#define DataFormats_CaloRecHit_TICLLayerClustersDeviceCollection_h

#include "DataFormats/CaloRecHit/interface/TICLLayerClusters.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using TICLLayerClustersDeviceCollection = PortableCollection<reco::TICLLayerClusters>;
}

#endif
