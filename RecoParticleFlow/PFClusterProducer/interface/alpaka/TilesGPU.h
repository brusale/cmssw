#ifndef LayerTilesGPU_h
#define LayerTilesGPU_h

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
// GPU Add
#include <cuda.h>
#include <cuda_runtime.h>

#include "GPUVecArray.h"
#include "TilesConstants.h"

// The type T is used to pass the number of bins in each dimension and the
// allowed ranges spanned. Anchillary quantitied, like the inverse of the bin
// width should also be provided. Code will not compile if any such information
// is missing.
template <typename T>
class TilesGPU {
 public:
  // constructor
  TilesGPU(){};

  __device__ void fill(float x, float y, int i) {
    tiles_[getGlobalBin(x, y)].push_back(i);
  }

  __host__ __device__ int getDim1Bin(float x) const {
    int dim1Bin = (x - T::minDim1) * T::invDim1BinSize;
    dim1Bin = std::clamp(dim1Bin, 0, T::nColumns - 1);
    return dim1Bin;
  }

  __host__ __device__ int getDim2Bin(float y) const {
    int dim2Bin = (y - T::minDim2) * T::invDim2BinSize;
    dim2Bin = std::clamp(dim2Bin, 0, T::nRows - 1);
    return dim2Bin;
  }

  __host__ __device__ int getGlobalBin(float x, float y) const {
    return getDim1Bin(x) + getDim2Bin(y) * T::nColumns;
  }

  __host__ __device__ int getGlobalBinByBin(int dim1_bin, int dim2_bin) const {
    return dim1_bin + dim2_bin * T::nColumns;
  }

  __host__ __device__ int4 searchBox(float dim1_min, float dim1_max,
                                     float dim2_min, float dim2_max) {
    int Bin1Min = getDim1Bin(dim1_min);
    int Bin1Max = getDim1Bin(dim1_max);
    int Bin2Min = getDim2Bin(dim2_min);
    int Bin2Max = getDim2Bin(dim2_max);
    return int4{Bin1Min, Bin1Max, Bin2Min, Bin2Max};
  }

  __host__ __device__ void clear() {
    for (auto& t : tiles_) t.reset();
  }

  __host__ __device__ GPU::VecArray<int, T::maxTileDepth>& operator[](
      int globalBinId) {
    return tiles_[globalBinId];
  }

 private:
  GPU::VecArray<GPU::VecArray<int, T::maxTileDepth>, T::nTiles> tiles_;
};
#endif
