#ifndef LayerTilesAlpaka_h
#define LayerTilesAlpaka_h

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>

#include "GPUVecArrayAlpaka.h"
#include "TilesConstants.h"

#if !defined(ALPAKA_ACC_GPU_CUDA_ENABLED) && \
    !defined(ALPAKA_ACC_GPU_HIP_ENABLED)
struct int4 {
  int x, y, z, w;
};
#endif

template <typename Acc, typename T>
class TilesAlpaka {
 public:
  using GPUVect = GPUAlpaka::VecArray<unsigned int, T::maxTileDepth>;
  // constructor
  TilesAlpaka(const Acc& acc) { acc_ = acc; };

  ALPAKA_FN_ACC
  void fill(float x, float y, int i) {
    tiles_[getGlobalBin(x, y)].push_back(acc_, i);
  }

  ALPAKA_FN_HOST_ACC int getDim1Bin(float x) const {
    int dim1Bin = (x - T::minDim1) * T::invDim1BinSize;
    // Equivalent of std::clamp(dim1Bin, 0, T::nColumns - 1);
    // which doesn't compile with gcc14 due to reference to __glibcxx_assert
    // See https://github.com/llvm/llvm-project/issues/95183
    // and https://gcc.gnu.org/bugzilla/show_bug.cgi?id=115740
    int new_max = std::max(dim1Bin, 0);
    dim1Bin = std::min(T::nColumns - 1, new_max);
    return dim1Bin;
  }

  ALPAKA_FN_HOST_ACC int getDim2Bin(float y) const {
    int dim2Bin = (y - T::minDim2) * T::invDim2BinSize;
    // See comment in getDim1Bin
    int new_max = std::max(dim2Bin, 0);
    dim2Bin = std::min(T::nRows - 1, new_max);
    return dim2Bin;
  }

  ALPAKA_FN_HOST_ACC int getDim2BinPhi(float phi) const {
    constexpr float r = T::nRows * M_1_PI * 0.5f;
    if (std::abs(phi) > M_PI) {
      float n = std::round(0.5f * phi / M_PI);
      phi *= (1 - 0.5f * n / M_PI);
    }
    int dim2Bin = (phi + M_PI) * r;
    return dim2Bin;
  }

  ALPAKA_FN_HOST_ACC int getGlobalBin(float x, float y) const {
    if (!T::phi) 
      return getDim1Bin(x) + getDim2Bin(y) * T::nColumns;
    else
      return getDim1Bin(x) + getDim2BinPhi(y) * T::nColumns;
  }

  ALPAKA_FN_HOST_ACC int getGlobalBinByBin(int dim1_bin, int dim2_bin) const {
    return dim1_bin + dim2_bin * T::nColumns;
  }

  ALPAKA_FN_HOST_ACC int4 searchBox(float dim1_min, float dim1_max,
                                    float dim2_min, float dim2_max) {
    int Bin1Min = getDim1Bin(dim1_min);
    int Bin1Max = getDim1Bin(dim1_max);
    int Bin2Min = getDim2Bin(dim2_min);
    int Bin2Max = getDim2Bin(dim2_max);
    return int4{Bin1Min, Bin1Max, Bin2Min, Bin2Max};
  }

  ALPAKA_FN_HOST_ACC void clear() {
    for (auto& t : tiles_) t.reset();
  }

  ALPAKA_FN_HOST_ACC void sort_unsafe(int i) {
    //for (int i = 0; i < T::nTiles; ++i)
      tiles_[i].sort_unsafe(acc_);
  }

  ALPAKA_FN_HOST_ACC GPUVect& operator[](int globalBinId) {
    return tiles_[globalBinId];
  }

 private:
  GPUAlpaka::VecArray<GPUVect, T::nTiles> tiles_;
  const Acc& acc_;
};
#endif
