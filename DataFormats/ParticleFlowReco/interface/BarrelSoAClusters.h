#ifndef DataFormats_ParticleFlowReco_interface_BarrelSoAClusters_h
#define DataFormats_ParticleFlowReco_interface_BarrelSoAClusters_h

#include <Eigen/Core>

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

GENERATE_SOA_LAYOUT(BarrelSoAClustersLayout,
                    // columns: one value per element
                    SOA_COLUMN(float, x),
                    SOA_COLUMN(float, y),
                    SOA_COLUMN(float, z),
                    SOA_COLUMN(float, energy),
                    SOA_COLUMN(int, cells),  // number of hits in the cluster
                    SOA_COLUMN(int, seed)    // This is the index of the seed of each cluster inside the RecHit SoA
)

using BarrelSoAClusters = BarrelSoAClustersLayout<>;

#endif
