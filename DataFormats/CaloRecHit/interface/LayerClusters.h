#ifndef DataFormats_CaloRecHit_LayerClusters_h
#define DataFormats_CaloRecHit_LayerClusters_h

// LayerClusters data format
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"
#include "DataFormats/SoATemplate/interface/SoACommon.h"

namespace reco {
    
  GENERATE_SOA_LAYOUT(LayerClustersSoA,
		      SOA_COLUMN(float, x),
		      SOA_COLUMN(float, y),
		      SOA_COLUMN(float, z), 
		      SOA_COLUMN(float, eta),
		      SOA_COLUMN(float, phi),
		      SOA_COLUMN(float, energy),
		      SOA_COLUMN(float, time),
		      SOA_COLUMN(float, error),
		      SOA_COLUMN(int, layerId), 
		      SOA_COLUMN(int, clusterIndex)
  )

  using LayerClusters = LayerClustersSoA<>;
  using LayerClustersView = LayerClusters::View;
  using LayerClustersConstView = LayerClusters::ConstView;

} //namespace reco

#endif
