#ifndef DataFormats_CaloRecHit_TICLLayerClustersHostCollection_h
#define DataFormats_CaloRecHit_TICLLayerClustersHostCollection_h

#include "DataFormats/CaloRecHit/interface/TICLLayerClusters.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"

namespace reco {
  using TICLLayerClustersHostCollection = PortableHostCollection<reco::TICLLayerClusters>;
}

#endif
