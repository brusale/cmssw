#ifndef DataFormats_HGCalReco_LayerClustersSoA_h
#define DataForamts_HGCalReco_LayerClustersSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

GENERATE_SOA_LAYOUT(LayerClustersSoALayout,
                    SOA_COLUMN(float, x),
                    SOA_COLUMN(float, y),
                    SOA_COLUMN(float, z),
                    SOA_COLUMN(float, eta),
                    SOA_COLUMN(float, phi),
                    SOA_COLUMN(float, energy),
                    SOA_COLUMN(float, error),
                    SOA_COLUMN(float, r_over_absz),
                    SOA_COLUMN(uint32_t, seed),
                    SOA_COLUMN(int, cells),
                    SOA_COLUMN(int, layerId),
                    SOA_COLUMN(int, clusterIndex),
                    SOA_COLUMN(int, algoId),
                    SOA_COLUMN(bool, isSilicon)
)

using LayerClustersSoA = LayerClustersSoALayout<>;
using LayerClustersSoAView = LayerClustersSoA::View;
using LayerClustersSoAConstView = LayerClustersSoA::ConstView;

#endif
