#ifndef DataFormats_ParticleFlowReco_interface_PFRecHitExtra_h
#define DataFormats_ParticleFlowReco_interface_PFRecHitExtra_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace reco {
  // SoA layout with delta, rho, weight, nearestHigher, clusterIndex, layer, isSeed, and cellsCount fields
  GENERATE_SOA_LAYOUT(PFRecHitExtraLayout,
                      // columns: one value per element
                      SOA_COLUMN(float, delta),
                      SOA_COLUMN(float, rho),
                      SOA_COLUMN(unsigned int, nearestHigher),
                      SOA_COLUMN(int, clusterIndex),
                      SOA_COLUMN(uint8_t, isSeed),
                      SOA_SCALAR(unsigned int, numberOfClustersScalar))
  
  using PFRecHitExtra = PFRecHitExtraLayout<>;
} // namespace reco
#endif
