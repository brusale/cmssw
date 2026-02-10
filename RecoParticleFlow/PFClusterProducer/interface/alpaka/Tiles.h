#ifndef LayerTiles_h
#define LayerTiles_h

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "TilesConstants.h"

// The type T is used to pass the number of bins in each dimension and the
// allowed ranges spanned. Anchillary quantitied, like the inverse of the bin
// width should also be provided. Code will not compile if any such information
// is missing.
template <typename T>
class Tiles {
 public:
  typedef T type;

  Tiles() { tiles_.resize(T::nTiles); }

  void fill(const std::vector<float>& x, const std::vector<float>& y) {
    auto cellsSize = x.size();
    for (unsigned int i = 0; i < cellsSize; ++i) {
      tiles_[getGlobalBin(x[i], y[i])].push_back(i);
    }
  }

  void fill(float x, float y, int i) {
    tiles_[getGlobalBin(x, y)].push_back(i);
  }

  int getDim1Bin(float x) const {
    int dim1Bin = (x - T::minDim1) * T::invDim1BinSize;
    dim1Bin = std::clamp(dim1Bin, 0, T::nColumns - 1);
    return dim1Bin;
  }

  int getDim2Bin(float y) const {
    int dim2Bin = (y - T::minDim2) * T::invDim2BinSize;
    dim2Bin = std::clamp(dim2Bin, 0, T::nRows - 1);
    return dim2Bin;
  }

  int getGlobalBin(float x, float y) const {
    return getDim1Bin(x) + getDim2Bin(y) * T::nColumns;
  }

  int getGlobalBinByBin(int dim1_bin, int dim2_bin) const {
    return dim1_bin + dim2_bin * T::nColumns;
  }

  std::array<int, 4> searchBox(float dim1_min, float dim1_max, float dim2_min,
                               float dim2_max) {
    int Bin1Min = getDim1Bin(dim1_min);
    int Bin1Max = getDim1Bin(dim1_max);
    int Bin2Min = getDim2Bin(dim2_min);
    int Bin2Max = getDim2Bin(dim2_max);
    return std::array<int, 4>({{Bin1Min, Bin1Max, Bin2Min, Bin2Max}});
  }

  void clear() {
    for (auto& t : tiles_) {
      t.clear();
    }
  }

  std::vector<int>& operator[](int globalBinId) { return tiles_[globalBinId]; }

 private:
  std::vector<std::vector<int>> tiles_;
};

#endif  // LayerTiles_h
