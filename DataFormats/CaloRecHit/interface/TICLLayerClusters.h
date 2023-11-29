#ifndef DataFormats_CaloRecHit_TICLLayerClusters_h
#define DataFormats_CaloRecHit_TICLLayerClusters_h

// LayerClusters data format
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"
#include "DataFormats/SoATemplate/interface/SoACommon.h"

#include "HeterogeneousCore/AlpakaInterface/interface/OneToManyAssoc.h"

struct HitAndFraction {
  float fraction;
  uint32_t hit;
};

using OneToManyHitAndFraction = cms::alpakatools::OneToManyAssoc<HitAndFraction, -1, -1>;

namespace reco {
  
  GENERATE_SOA_LAYOUT(TICLLayerClustersSoA,
		      SOA_COLUMN(float, x),
		      SOA_COLUMN(float, y),
		      SOA_COLUMN(float, z), 
		      SOA_COLUMN(float, eta),
		      SOA_COLUMN(float, phi),
		      SOA_COLUMN(float, time),
		      SOA_COLUMN(float, err_pos),
          SOA_COLUMN(float, energy),
		      SOA_SCALAR(OneToManyHitAndFraction, hitsAndFractions)
  )

  using TICLLayerClusters = TICLLayerClustersSoA<>;
  using TICLLayerClustersView = TICLLayerClusters::View;
  using TICLLayerClustersConstView = TICLLayerClusters::ConstView;

} //namespace reco

#endif
