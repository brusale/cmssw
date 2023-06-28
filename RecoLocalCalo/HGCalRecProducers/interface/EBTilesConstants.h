// Authors: Felice Pantaleo - felice.pantaleo@cern.ch
// Date: 03/2019

#ifndef RecoLocalCalo_HGCalRecProducer_EBTilesConstants_h
#define RecoLocalCalo_HGCalRecProducer_EBTilesConstants_h

#include "DataFormats/Math/interface/constexpr_cmath.h"
#include <cmath>
#include <cstdint>
#include <array>

struct EBTilesConstants {
  static constexpr float tileSize = 5.f;
  static constexpr float minX = -285.f;
  static constexpr float maxX = 285.f;
  static constexpr float minY = -285.f;
  static constexpr float maxY = 285.f;
  static constexpr int nColumns = reco::ceil((maxX - minX) / tileSize);
  static constexpr int nRows = reco::ceil((maxY - minY) / tileSize);
  static constexpr float tileSizeEtaPhi =3*0.0175f;
  static constexpr float minEta = -1.5f;
  static constexpr float maxEta = 1.5f;
  static constexpr int nColumnsEta = reco::ceil((maxEta - minEta) / tileSizeEtaPhi);
  static constexpr int nRowsPhi = reco::ceil(2. * M_PI / tileSizeEtaPhi);
  static constexpr int nTiles = nColumns * nRows + nColumnsEta * nRowsPhi;
  static constexpr float cellWidthEta = 0.0175f;
  static constexpr float cellWidthPhi = 0.0175f;
};

#endif
