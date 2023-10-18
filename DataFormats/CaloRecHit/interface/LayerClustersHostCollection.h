#ifndef DataFormats_CaloRecHit_LayerClustersHostCollection_h
#define DataFormats_CaloRecHit_LayerClustersHostCollection_h

#include "DataFormats/CaloRecHit/interface/LayerClusters.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"

namespace reco {
  using LayerClustersHostCollection = PortableHostCollection<reco::LayerClusters>;
}

#endif
