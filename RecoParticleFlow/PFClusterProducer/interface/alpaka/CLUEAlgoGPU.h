#ifndef CLUEAlgoGPU_h
#define CLUEAlgoGPU_h
#include <math.h>

#include <iostream>
#include <limits>

// GPU Add
#include <cuda.h>
#include <cuda_runtime.h>
// for timing
#include <chrono>
#include <ctime>

#include "CLUEAlgo.h"
#include "TilesGPU.h"

#define CHECK_CUDA_ERROR(val) check((val), #val, __FILE__, __LINE__)
template <typename T>
void check(T err, const char *const func, const char *const file,
           const int line) {
  if (err != cudaSuccess) {
    std::cerr << "CUDA Runtime Error at: " << file << ":" << line << std::endl;
    std::cerr << cudaGetErrorString(err) << " " << func << std::endl;

    std::exit(EXIT_FAILURE);
  }
}

class WorkDivByPoints;
class WorkDivByTile;

static const int maxNSeeds = 100000;
static const int maxNFollowers = 32;
static const int localStackSizePerSeed = 32;

struct PointsPtr {
  float *x;
  float *y;
  int *layer;
  float *weight;
  float *sigmaNoise;

  float *rho;
  float *delta;
  int *nearestHigher;
  int *clusterIndex;
  uint8_t *isSeed;
};

// The type T is used to pass the number of bins in each dimension and the
// allowed ranges spanned. Anchillary quantitied, like the inverse of the bin
// width should also be provided. Code will not compile if any such information
// is missing.
template <typename T, int NLAYERS, typename W=WorkDivByPoints>
class CLUEAlgoGPU : public CLUEAlgo<T, NLAYERS> {
  // inherit from CLUEAlgo

 public:
  // constructor
  CLUEAlgoGPU(float dc, float rhoc, float outlierDeltaFactor, bool verbose, bool useAbsoluteSigma=false)
      : CLUEAlgo<T, NLAYERS>(dc, rhoc, outlierDeltaFactor, verbose, useAbsoluteSigma) {
    init_device();
  }
  // destructor
  ~CLUEAlgoGPU() { free_device(); }

  // public methods
  void makeClusters();  // overwrite base class

  // Bring base-class public variables into the scope of this template derived
  // class
  using CLUEAlgo<T, NLAYERS>::dc_;
  using CLUEAlgo<T, NLAYERS>::rhoc_;
  using CLUEAlgo<T, NLAYERS>::kappa_;
  using CLUEAlgo<T, NLAYERS>::outlierDeltaFactor_;
  using CLUEAlgo<T, NLAYERS>::verbose_;
  using CLUEAlgo<T, NLAYERS>::points_;
  using CLUEAlgo<T, NLAYERS>::useAbsoluteSigma_;

 private:
  // private variables

  cudaStream_t stream_;
  // algorithm internal variables
  PointsPtr d_points;
  TilesGPU<T> *d_hist;
  GPU::VecArray<int, maxNSeeds> *d_seeds;
  GPU::VecArray<int, maxNFollowers> *d_followers;

  // private methods
  void init_device() {
    // Create our own cuda stream
    CHECK_CUDA_ERROR(cudaStreamCreate(&stream_));
    // Allocate memory
    unsigned int reserve = 1000000;
    // input variables
    CHECK_CUDA_ERROR(
        cudaMallocAsync(&d_points.x, sizeof(float) * reserve, stream_));
    CHECK_CUDA_ERROR(
        cudaMallocAsync(&d_points.y, sizeof(float) * reserve, stream_));
    CHECK_CUDA_ERROR(
        cudaMallocAsync(&d_points.layer, sizeof(int) * reserve, stream_));
    CHECK_CUDA_ERROR(
        cudaMallocAsync(&d_points.weight, sizeof(float) * reserve, stream_));
    if (useAbsoluteSigma_)
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&d_points.sigmaNoise, sizeof(float) * reserve, stream_));
    // result variables
    CHECK_CUDA_ERROR(
        cudaMallocAsync(&d_points.rho, sizeof(float) * reserve, stream_));
    CHECK_CUDA_ERROR(
        cudaMallocAsync(&d_points.delta, sizeof(float) * reserve, stream_));
    CHECK_CUDA_ERROR(cudaMallocAsync(&d_points.nearestHigher,
                                     sizeof(int) * reserve, stream_));
    CHECK_CUDA_ERROR(cudaMallocAsync(&d_points.clusterIndex,
                                     sizeof(int) * reserve, stream_));
    CHECK_CUDA_ERROR(
        cudaMallocAsync(&d_points.isSeed, sizeof(uint8_t) * reserve, stream_));
    // algorithm internal variables
    CHECK_CUDA_ERROR(
        cudaMallocAsync(&d_hist, sizeof(TilesGPU<T>) * NLAYERS, stream_));
    CHECK_CUDA_ERROR(cudaMallocAsync(
        &d_seeds, sizeof(GPU::VecArray<int, maxNSeeds>), stream_));
    CHECK_CUDA_ERROR(cudaMallocAsync(
        &d_followers, sizeof(GPU::VecArray<int, maxNFollowers>) * reserve,
        stream_));
  }

  void free_device() {
    // input variables
    CHECK_CUDA_ERROR(cudaFreeAsync(d_points.x, stream_));
    CHECK_CUDA_ERROR(cudaFreeAsync(d_points.y, stream_));
    CHECK_CUDA_ERROR(cudaFreeAsync(d_points.layer, stream_));
    CHECK_CUDA_ERROR(cudaFreeAsync(d_points.weight, stream_));
    if (useAbsoluteSigma_)
      CHECK_CUDA_ERROR(cudaFreeAsync(d_points.sigmaNoise, stream_));
    // result variables
    CHECK_CUDA_ERROR(cudaFreeAsync(d_points.rho, stream_));
    CHECK_CUDA_ERROR(cudaFreeAsync(d_points.delta, stream_));
    CHECK_CUDA_ERROR(cudaFreeAsync(d_points.nearestHigher, stream_));
    CHECK_CUDA_ERROR(cudaFreeAsync(d_points.clusterIndex, stream_));
    CHECK_CUDA_ERROR(cudaFreeAsync(d_points.isSeed, stream_));
    // algorithm internal variables
    CHECK_CUDA_ERROR(cudaFreeAsync(d_hist, stream_));
    CHECK_CUDA_ERROR(cudaFreeAsync(d_seeds, stream_));
    CHECK_CUDA_ERROR(cudaFreeAsync(d_followers, stream_));

    CHECK_CUDA_ERROR(cudaStreamDestroy(stream_));
  }

  void copy_todevice() {
    // input variables
    CHECK_CUDA_ERROR(cudaMemcpyAsync(d_points.x, points_.p_x,
                                     sizeof(float) * points_.n,
                                     cudaMemcpyHostToDevice, stream_));
    CHECK_CUDA_ERROR(cudaMemcpyAsync(d_points.y, points_.p_y,
                                     sizeof(float) * points_.n,
                                     cudaMemcpyHostToDevice, stream_));
    CHECK_CUDA_ERROR(cudaMemcpyAsync(d_points.layer, points_.p_layer,
                                     sizeof(int) * points_.n,
                                     cudaMemcpyHostToDevice, stream_));
    CHECK_CUDA_ERROR(cudaMemcpyAsync(d_points.weight, points_.p_weight,
                                     sizeof(float) * points_.n,
                                     cudaMemcpyHostToDevice, stream_));
    if (useAbsoluteSigma_)
      CHECK_CUDA_ERROR(cudaMemcpyAsync(d_points.sigmaNoise, points_.p_sigmaNoise,
                                       sizeof(float) * points_.n,
                                       cudaMemcpyHostToDevice, stream_));
  }

  void clear_internal_buffers() {
    // // result variables
    CHECK_CUDA_ERROR(cudaMemsetAsync(d_points.rho, 0x00,
                                     sizeof(float) * points_.n, stream_));
    CHECK_CUDA_ERROR(cudaMemsetAsync(d_points.delta, 0x00,
                                     sizeof(float) * points_.n, stream_));
    CHECK_CUDA_ERROR(cudaMemsetAsync(d_points.nearestHigher, 0x00,
                                     sizeof(int) * points_.n, stream_));
    CHECK_CUDA_ERROR(cudaMemsetAsync(d_points.clusterIndex, 0x00,
                                     sizeof(int) * points_.n, stream_));
    CHECK_CUDA_ERROR(cudaMemsetAsync(d_points.isSeed, 0x00,
                                     sizeof(uint8_t) * points_.n, stream_));
    // algorithm internal variables
    CHECK_CUDA_ERROR(
        cudaMemsetAsync(d_hist, 0x00, sizeof(TilesGPU<T>) * NLAYERS, stream_));
    CHECK_CUDA_ERROR(cudaMemsetAsync(
        d_seeds, 0x00, sizeof(GPU::VecArray<int, maxNSeeds>), stream_));
    CHECK_CUDA_ERROR(cudaMemsetAsync(
        d_followers, 0x00,
        sizeof(GPU::VecArray<int, maxNFollowers>) * points_.n, stream_));
  }

  void copy_tohost() {
    // result variables
    CHECK_CUDA_ERROR(cudaMemcpyAsync(
        points_.clusterIndex.data(), d_points.clusterIndex,
        sizeof(int) * points_.n, cudaMemcpyDeviceToHost, stream_));
    if (verbose_) {
      // other variables, copy only when verbose_==True
      CHECK_CUDA_ERROR(cudaMemcpyAsync(points_.rho.data(), d_points.rho,
                                       sizeof(float) * points_.n,
                                       cudaMemcpyDeviceToHost, stream_));
      CHECK_CUDA_ERROR(cudaMemcpyAsync(points_.delta.data(), d_points.delta,
                                       sizeof(float) * points_.n,
                                       cudaMemcpyDeviceToHost, stream_));
      CHECK_CUDA_ERROR(cudaMemcpyAsync(
          points_.nearestHigher.data(), d_points.nearestHigher,
          sizeof(int) * points_.n, cudaMemcpyDeviceToHost, stream_));
      CHECK_CUDA_ERROR(cudaMemcpyAsync(points_.isSeed.data(), d_points.isSeed,
                                       sizeof(uint8_t) * points_.n,
                                       cudaMemcpyDeviceToHost, stream_));
    }
  }

  // #endif // __CUDACC__
};

template <typename T>
__global__ void kernel_compute_histogram(TilesGPU<T> *d_hist,
                                         const PointsPtr d_points,
                                         int numberOfPoints) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < numberOfPoints) {
    // push index of points into tiles
    d_hist[d_points.layer[i]].fill(d_points.x[i], d_points.y[i], i);
  }
}  // kernel kernel_compute_histogram

template <typename T>
__global__ void kernel_calculate_densityTile(TilesGPU<T> *d_hist,
                                         PointsPtr d_points, float dc,
                                         int numberOfPoints) {
  int layeri = blockIdx.y;
  int globalBinOnLayer = blockIdx.x;
  int bin = threadIdx.x;

  if (bin < d_hist[layeri][globalBinOnLayer].size()) {
    int i = d_hist[layeri][globalBinOnLayer][bin];
    float rhoi{0.};
    int layeri = d_points.layer[i];
    float xi = d_points.x[i];
    float yi = d_points.y[i];

    // get search box
    int4 search_box =
        d_hist[layeri].searchBox(xi - dc, xi + dc, yi - dc, yi + dc);

    // loop over bins in the search box
    for (int xBin = search_box.x; xBin < search_box.y + 1; ++xBin) {
      for (int yBin = search_box.z; yBin < search_box.w + 1; ++yBin) {
        // get the id of this bin
        int binId = d_hist[layeri].getGlobalBinByBin(xBin, yBin);
        // get the size of this bin
        int binSize = d_hist[layeri][binId].size();

        // interate inside this bin
        for (int binIter = 0; binIter < binSize; binIter++) {
          int j = d_hist[layeri][binId][binIter];
          // query N_{dc_}(i)
          float xj = d_points.x[j];
          float yj = d_points.y[j];
          float dist_ij_2 =
              ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj));
          rhoi += (dist_ij_2 <= dc*dc) * (i == j ? 1.f : 0.5f) * d_points.weight[j];
          /*
          if (dist_ij_2 <= dc*dc) {
            // sum weights within N_{dc_}(i)
            rhoi += (i == j ? 1.f : 0.5f) * d_points.weight[j];
          }
          */
        }  // end of interate inside this bin
      }
    }  // end of loop over bins in search box
    d_points.rho[i] = rhoi;
  }
}  // kernel kernel_calculate_density

template <typename T>
__global__ void kernel_calculate_density(TilesGPU<T> *d_hist,
                                         PointsPtr d_points, float dc,
                                         int numberOfPoints) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < numberOfPoints) {
    float rhoi{0.};
    int layeri = d_points.layer[i];
    float xi = d_points.x[i];
    float yi = d_points.y[i];

    // get search box
    int4 search_box =
        d_hist[layeri].searchBox(xi - dc, xi + dc, yi - dc, yi + dc);

    // loop over bins in the search box
    for (int xBin = search_box.x; xBin < search_box.y + 1; ++xBin) {
      for (int yBin = search_box.z; yBin < search_box.w + 1; ++yBin) {
        // get the id of this bin
        int binId = d_hist[layeri].getGlobalBinByBin(xBin, yBin);
        // get the size of this bin
        int binSize = d_hist[layeri][binId].size();

        // interate inside this bin
        for (int binIter = 0; binIter < binSize; binIter++) {
          int j = d_hist[layeri][binId][binIter];
          // query N_{dc_}(i)
          float xj = d_points.x[j];
          float yj = d_points.y[j];
          float dist_ij_2 =
              ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj));
          rhoi += (dist_ij_2 <= dc*dc) * (i == j ? 1.f : 0.5f) * d_points.weight[j];
          /*
          if (dist_ij_2 <= dc*dc) {
            // sum weights within N_{dc_}(i)
            rhoi += (i == j ? 1.f : 0.5f) * d_points.weight[j];
          }
          */
        }  // end of interate inside this bin
      }
    }  // end of loop over bins in search box
    d_points.rho[i] = rhoi;
  }
}  // kernel kernel_calculate_density

template <typename T>
__global__ void kernel_calculate_distanceToHigher(TilesGPU<T> *d_hist,
                                                  PointsPtr d_points,
                                                  float outlierDeltaFactor,
                                                  float dc, int numberOfPoints) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  const float dm = outlierDeltaFactor * dc;

  if (i < numberOfPoints) {
    int layeri = d_points.layer[i];

    float deltai = std::numeric_limits<float>::max();
    int nearestHigheri = -1;
    float xi = d_points.x[i];
    float yi = d_points.y[i];
    float rhoi = d_points.rho[i];

    // get search box
    int4 search_box =
        d_hist[layeri].searchBox(xi - dm, xi + dm, yi - dm, yi + dm);

    // loop over all bins in the search box
    for (int xBin = search_box.x; xBin < search_box.y + 1; ++xBin) {
      for (int yBin = search_box.z; yBin < search_box.w + 1; ++yBin) {
        // get the id of this bin
        int binId = d_hist[layeri].getGlobalBinByBin(xBin, yBin);
        // get the size of this bin
        int binSize = d_hist[layeri][binId].size();

        // interate inside this bin
        for (int binIter = 0; binIter < binSize; binIter++) {
          int j = d_hist[layeri][binId][binIter];
          // query N'_{dm}(i)
          float xj = d_points.x[j];
          float yj = d_points.y[j];
          float dist_ij_2 =
              ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj));
          bool foundHigher = (d_points.rho[j] > rhoi);
          // in the rare case where rho is the same, use detid
          foundHigher = foundHigher || ((d_points.rho[j] == rhoi) && (j > i));
          if (foundHigher && dist_ij_2 <= dm*dm) {  // definition of N'_{dm}(i)
            // find the nearest point within N'_{dm}(i)
            if (dist_ij_2 < deltai) {
              // update deltai and nearestHigheri
              deltai = dist_ij_2;
              nearestHigheri = j;
            }
          }
        }  // end of interate inside this bin
      }
    }  // end of loop over bins in search box
    d_points.delta[i] = std::sqrt(deltai);
    d_points.nearestHigher[i] = nearestHigheri;
  }
}  // kernel kernel_calculate_distanceToHigher

template <typename T>
__global__ void kernel_calculate_distanceToHigherTile(TilesGPU<T> *d_hist,
                                                  PointsPtr d_points,
                                                  float outlierDeltaFactor,
                                                  float dc, int numberOfPoints) {
  int layeri = blockIdx.y;
  int globalBinOnLayer = blockIdx.x;
  int bin = threadIdx.x;

  if (bin < d_hist[layeri][globalBinOnLayer].size()) {
    const float dm = outlierDeltaFactor * dc;

    int i = d_hist[layeri][globalBinOnLayer][bin];
    float deltai = std::numeric_limits<float>::max();
    int nearestHigheri = -1;
    float xi = d_points.x[i];
    float yi = d_points.y[i];
    float rhoi = d_points.rho[i];

    // get search box
    int4 search_box =
        d_hist[layeri].searchBox(xi - dm, xi + dm, yi - dm, yi + dm);

    // loop over all bins in the search box
    for (int xBin = search_box.x; xBin < search_box.y + 1; ++xBin) {
      for (int yBin = search_box.z; yBin < search_box.w + 1; ++yBin) {
        // get the id of this bin
        int binId = d_hist[layeri].getGlobalBinByBin(xBin, yBin);
        // get the size of this bin
        int binSize = d_hist[layeri][binId].size();

        // interate inside this bin
        for (int binIter = 0; binIter < binSize; binIter++) {
          int j = d_hist[layeri][binId][binIter];
          // query N'_{dm}(i)
          float xj = d_points.x[j];
          float yj = d_points.y[j];
          float dist_ij_2 =
              ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj));
          bool foundHigher = (d_points.rho[j] > rhoi);
          // in the rare case where rho is the same, use detid
          foundHigher = foundHigher || ((d_points.rho[j] == rhoi) && (j > i));
          if (foundHigher && dist_ij_2 <= dm*dm) {  // definition of N'_{dm}(i)
            // find the nearest point within N'_{dm}(i)
            if (dist_ij_2 < deltai) {
              // update deltai and nearestHigheri
              deltai = dist_ij_2;
              nearestHigheri = j;
            }
          }
        }  // end of interate inside this bin
      }
    }  // end of loop over bins in search box
    d_points.delta[i] = std::sqrt(deltai);
    d_points.nearestHigher[i] = nearestHigheri;
  }
}  // kernel kernel_calculate_distanceToHigher

__global__ void kernel_find_clusters(
    GPU::VecArray<int, maxNSeeds> *d_seeds,
    GPU::VecArray<int, maxNFollowers> *d_followers, PointsPtr d_points,
    float outlierDeltaFactor, float dc, float rhoc, int numberOfPoints) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < numberOfPoints) {
    // initialize clusterIndex
    d_points.clusterIndex[i] = -1;
    // determine seed or outlier
    float deltai = d_points.delta[i];
    float rhoi = d_points.rho[i];
    bool isSeed = (deltai > dc) && (rhoi >= rhoc);
    bool isOutlier = (deltai > outlierDeltaFactor * dc) && (rhoi < rhoc);

    if (isSeed) {
      // set isSeed as 1
      d_points.isSeed[i] = 1;
      d_seeds[0].push_back(i);  // head of d_seeds
    } else {
      if (!isOutlier) {
        assert(d_points.nearestHigher[i] < numberOfPoints);
        // register as follower at its nearest higher
        d_followers[d_points.nearestHigher[i]].push_back(i);
      }
    }
  }
}  // kernel kernel_find_clusters

__global__ void kernel_find_clusters_kappa(
    GPU::VecArray<int, maxNSeeds> *d_seeds,
    GPU::VecArray<int, maxNFollowers> *d_followers, PointsPtr d_points,
    float outlierDeltaFactor, float dc, float kappa, int numberOfPoints) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < numberOfPoints) {
    // initialize clusterIndex
    d_points.clusterIndex[i] = -1;
    // determine seed or outlier
    float deltai = d_points.delta[i];
    float rhoi = d_points.rho[i];
    float rhoc = d_points.sigmaNoise[i] * kappa;
    bool isSeed = (deltai > dc) && (rhoi >= rhoc);
    bool isOutlier = (deltai > outlierDeltaFactor * dc) && (rhoi < rhoc);

    if (isSeed) {
      // set isSeed as 1
      d_points.isSeed[i] = 1;
      d_seeds[0].push_back(i);  // head of d_seeds
    } else {
      if (!isOutlier) {
        assert(d_points.nearestHigher[i] < numberOfPoints);
        // register as follower at its nearest higher
        d_followers[d_points.nearestHigher[i]].push_back(i);
      }
    }
  }
}  // kernel kernel_find_clusters

__global__ void kernel_assign_clusters(
    const GPU::VecArray<int, maxNSeeds> *d_seeds,
    const GPU::VecArray<int, maxNFollowers> *d_followers, PointsPtr d_points,
    int numberOfPoints) {
  int idxCls = blockIdx.x * blockDim.x + threadIdx.x;
  const auto &seeds = d_seeds[0];
  const auto nSeeds = seeds.size();
  if (idxCls < nSeeds) {
    int localStack[localStackSizePerSeed] = {-1};
    int localStackSize = 0;

    // asgine cluster to seed[idxCls]
    int idxThisSeed = seeds[idxCls];
    d_points.clusterIndex[idxThisSeed] = idxCls;
    // push_back idThisSeed to localStack
    localStack[localStackSize] = idxThisSeed;
    localStackSize++;
    // process all elements in localStack
    while (localStackSize > 0) {
      // get last element of localStack
      int idxEndOflocalStack = localStack[localStackSize - 1];

      int temp_clusterIndex = d_points.clusterIndex[idxEndOflocalStack];
      // pop_back last element of localStack
      localStack[localStackSize - 1] = -1;
      localStackSize--;

      // loop over followers of last element of localStack
      for (int j : d_followers[idxEndOflocalStack]) {
        // // pass id to follower
        d_points.clusterIndex[j] = temp_clusterIndex;
        // push_back follower to localStack
        localStack[localStackSize] = j;
        localStackSize++;
      }
    }
  }
}  // kernel kernel_assign_clusters

template <typename T, int NLAYERS, typename W>
void CLUEAlgoGPU<T, NLAYERS, W>::makeClusters() {
  copy_todevice();
  clear_internal_buffers();

  ////////////////////////////////////////////
  // calculate rho, delta and find seeds
  // 1 point per thread
  ////////////////////////////////////////////
  const dim3 blockSize(64, 1, 1);
  const dim3 gridSize(ceil(points_.n / static_cast<float>(blockSize.x)), 1, 1);
  kernel_compute_histogram<T>
      <<<gridSize, blockSize, 0, stream_>>>(d_hist, d_points, points_.n);

  if constexpr (std::is_same_v<W, WorkDivByPoints>) {
    kernel_calculate_density<T>
      <<<gridSize, blockSize, 0, stream_>>>(d_hist, d_points, dc_, points_.n);
    kernel_calculate_distanceToHigher<T><<<gridSize, blockSize, 0, stream_>>>(
        d_hist, d_points, outlierDeltaFactor_, dc_, points_.n);
  }

  if constexpr (std::is_same_v<W, WorkDivByTile>) {
    const dim3 gridSizeTile(T::nTiles, NLAYERS, 1);
    const dim3 blockSizeTile(T::maxTileDepth, 1, 1);
    kernel_calculate_densityTile<T>
      <<<gridSizeTile, blockSizeTile, 0, stream_>>>(d_hist, d_points, dc_, points_.n);
    kernel_calculate_distanceToHigherTile<T><<<gridSizeTile, blockSizeTile, 0, stream_>>>(
        d_hist, d_points, outlierDeltaFactor_, dc_, points_.n);
  }

  if (!useAbsoluteSigma_) {
    kernel_find_clusters<<<gridSize, blockSize, 0, stream_>>>(d_seeds, d_followers, d_points,
                                                  outlierDeltaFactor_, dc_, rhoc_,
                                                  points_.n);
  } else {
    kernel_find_clusters_kappa<<<gridSize, blockSize, 0, stream_>>>(d_seeds, d_followers, d_points,
                                                        outlierDeltaFactor_, dc_, kappa_,
                                                        points_.n);
  }

  ////////////////////////////////////////////
  // assign clusters
  // 1 point per seeds
  ////////////////////////////////////////////
  const dim3 gridSize_nseeds(ceil(maxNSeeds / 1024.0), 1, 1);
  kernel_assign_clusters<<<gridSize_nseeds, blockSize, 0, stream_>>>(
      d_seeds, d_followers, d_points, points_.n);

  copy_tohost();
  CHECK_CUDA_ERROR(cudaStreamSynchronize(stream_));
}

#endif
