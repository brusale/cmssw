#ifndef DataFormats_ParticleFlowReco_interface_PFRecHitExtraHostCollection_h
#define DataFormats_ParticleFlowReco_interface_PFRecHitExtraHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitExtra.h"

namespace reco {
  // SoA with delta, rho, weight, nearestHigher, clusterIndex, layer, isSeed, and cellsCount fields in host memory
  using PFRecHitExtraHostCollection = PortableHostCollection<reco::PFRecHitExtra>;
} // namespace reco

#endif  // DataFormats_ParticleFlowReco_interface_PFRecHitExtraHostCollection_h
