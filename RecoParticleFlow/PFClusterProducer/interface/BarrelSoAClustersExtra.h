#ifndef DataFormats_ParticleFlowReco_interface_BarrelSoAClustersExtra_h
#define DataFormats_ParticleFlowReco_interface_BarrelSoAClustersExtra_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

GENERATE_SOA_LAYOUT(
    BarrelSoAClustersExtraLayout,
    // columns: one value per element
    SOA_COLUMN(float, total_weight),
    SOA_COLUMN(float, total_weight_log),
    SOA_COLUMN(float, maxEnergyValue),
    SOA_COLUMN(int, maxEnergyIndex)  // Index in the RecHitSoA of the rechit with highest energy in each cluster
)

using BarrelSoAClustersExtra = BarrelSoAClustersExtraLayout<>;

#endif
