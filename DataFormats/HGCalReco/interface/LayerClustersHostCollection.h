#ifndef DataFormats_HGCalReco_LayerClustersHostCollection_h
#define DataFormats_HGCalReco_LayerClustersHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/HGCalReco/interface/LayerClustersSoA.h"

using LayerClustersCollection = PortableHostCollection<LayerClustersSoA>;
using LayerClustersCollectionView = PortableHostCollection<LayerClustersSoA>::View;
using LayerClustersCollectionConstView = PortableHostCollection<LayerClustersSoA>::ConstView;

#endif
