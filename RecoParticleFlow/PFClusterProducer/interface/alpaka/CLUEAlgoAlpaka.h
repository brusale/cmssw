#pragma once
#include <alpaka/alpaka.hpp>
#include <cfloat>
#include <chrono>
#include <optional>

#include "CLUEAlgo.h"
#include "TilesAlpaka.h"

#define ORDER_TILE 0

#define DECLARE_TASKTYPE_AND_KERNEL(ACC, NAME, ...)                            \
  struct Kernel##NAME {};                                                      \
  ALPAKA_FN_ACC void operator()(ACC const &acc, Kernel##NAME dummy,            \
                                ##__VA_ARGS__) const

// Main interface to select the proper accelerator **at compile time**.
namespace alpaka {
//! Alias for the default accelerator used by examples. From a list of
//! all accelerators the first one which is enabled is chosen.
//! AccCpuSerial is selected last.
template <class TDim, class TIdx>
#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
using SelectedAcc = alpaka::AccGpuCudaRt<TDim, TIdx>;
#elif defined(ALPAKA_ACC_GPU_HIP_ENABLED)
using SelectedAcc = alpaka::AccGpuHipRt<TDim, TIdx>;
#elif defined(ALPAKA_ACC_CPU_B_OMP2_T_SEQ_ENABLED)
using SelectedAcc = alpaka::AccCpuOmp2Blocks<TDim, TIdx>;
#elif defined(ALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED)
using SelectedAcc = alpaka::AccCpuTbbBlocks<TDim, TIdx>;
#elif defined(ALPAKA_ACC_CPU_B_SEQ_T_FIBERS_ENABLED)
using SelectedAcc = alpaka::AccCpuFibers<TDim, TIdx>;
#elif defined(ALPAKA_ACC_CPU_B_SEQ_T_OMP2_ENABLED)
using SelectedAcc = alpaka::AccCpuOmp2Threads<TDim, TIdx>;
#elif defined(ALPAKA_ACC_CPU_B_SEQ_T_THREADS_ENABLED)
using SelectedAcc = alpaka::AccCpuThreads<TDim, TIdx>;
#elif defined(ALPAKA_ACC_ANY_BT_OMP5_ENABLED)
using SelectedAcc = alpaka::AccOmp5<TDim, TIdx>;
#elif defined(ALPAKA_ACC_ANY_BT_OACC_ENABLED)
using SelectedAcc = alpaka::AccOacc<TDim, TIdx>;
#elif defined(ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED)
using SelectedAcc = alpaka::AccCpuSerial<TDim, TIdx>;
#else
class SelectedAcc;
#warning "No supported backend selected."
#endif
} // namespace alpaka

// Maximum number of uniques seeds that could be handled. A higher number of
// potential seed will trigger an exception.
static const int maxNSeeds = 262144;

// Maximum number of followers that could be handled. A higher number of
// followers will trigger an exception.
static const int maxNFollowers = 128;

// Maximum size of the local stack used to assign clusters to seeds and
// followers. It should be at least as big as the maximum allowed number of
// followers. Adding more elements with respect to the reserved size will
// trigger an exception.
static const int localStackSizePerSeed = 128;

// The type T is used to pass the number of bins in each dimension and the
// allowed ranges spanned. Anchillary quantitied, like the inverse of the bin
// width should also be provided. Code will not compile if any such information
// is missing.
template <typename TAcc, typename TQueue, typename T, int NLAYERS>
class CLUEAlgoAlpaka : public CLUEAlgo<T, NLAYERS> {
public:
  using Dim = alpaka::Dim<TAcc>;
  using Idx = alpaka::Idx<TAcc>;

  template <typename TT> using BufAccT = alpaka::Buf<TAcc, TT, Dim, Idx>;

  template <typename TT>
  using ViewHostT = alpaka::ViewPlainPtr<alpaka::DevCpu, TT, Dim, Idx>;

  using LayerTilesAcc = TilesAlpaka<TAcc, T>;
  using BufLayerTiles = BufAccT<LayerTilesAcc>;

  using BufVecArrSeeds = BufAccT<GPUAlpaka::VecArray<int, maxNSeeds>>;
  using BufVecArrFollowers = BufAccT<GPUAlpaka::VecArray<int, maxNFollowers>>;
  //
  // Bring base-class public variables into the scope of this template derived
  // class
  using CLUEAlgo<T, NLAYERS>::dc_;
  using CLUEAlgo<T, NLAYERS>::rhoc_;
  using CLUEAlgo<T, NLAYERS>::kappa_;
  using CLUEAlgo<T, NLAYERS>::outlierDeltaFactor_;
  using CLUEAlgo<T, NLAYERS>::verbose_;
  using CLUEAlgo<T, NLAYERS>::useAbsoluteSigma_;
  using CLUEAlgo<T, NLAYERS>::points_;

  struct PointsBuf {
    // Input Buffers
    std::optional<BufAccT<float>> x;
    std::optional<BufAccT<float>> y;
    std::optional<BufAccT<int>> layer;
    std::optional<BufAccT<float>> weight;
    std::optional<BufAccT<float>> sigmaNoise;

    // Output Buffers
    std::optional<BufAccT<float>> rho;
    std::optional<BufAccT<float>> delta;
    std::optional<BufAccT<unsigned int>> nearestHigher;
    std::optional<BufAccT<int>> clusterIndex;
    std::optional<BufAccT<uint8_t>> isSeed;
  };

  class DeviceRunner {
  public:
    struct DeviceRawPointers {
      // Pointers to Input buffers on device
      float *x;
      float *y;
      int *layer;
      float *weight;
      float *sigmaNoise;
      uint32_t *detid;

      // Pointers to output buffers on device
      float *rho;
      float *delta;
      unsigned int *nearestHigher;
      int *clusterIndex;
      uint8_t *isSeed;

      // LayerTiles and utility data structures
      TilesAlpaka<TAcc, T> *hist_;
      GPUAlpaka::VecArray<int, maxNSeeds> *seeds_;
      GPUAlpaka::VecArray<int, maxNFollowers> *followers_;
    };

    // Kernel and KernelTask definitions
    DECLARE_TASKTYPE_AND_KERNEL(TAcc, ComputeHistogram,
                                const unsigned int num_elements);
    DECLARE_TASKTYPE_AND_KERNEL(TAcc, SortHistogram);
    DECLARE_TASKTYPE_AND_KERNEL(TAcc, InitializeHistogram,
                                const unsigned int num_elements,
                                const float outlierDeltaFactor,
                                const float dc);
    DECLARE_TASKTYPE_AND_KERNEL(TAcc, ComputeLocalDensity, float dc,
                                const unsigned int num_elements);
    DECLARE_TASKTYPE_AND_KERNEL(TAcc, ComputeDistanceToHigherNoDetId,
                                float outlierDeltaFactor, float dc,
                                const unsigned int num_elements);
    DECLARE_TASKTYPE_AND_KERNEL(TAcc, ComputeDistanceToHigher,
                                float outlierDeltaFactor, float dc,
                                const unsigned int num_elements);
    DECLARE_TASKTYPE_AND_KERNEL(TAcc, FindClusters, float outlierDeltaFactor,
                                float dc, float rhoc,
                                const unsigned int num_elements);
    DECLARE_TASKTYPE_AND_KERNEL(TAcc, FindClustersKappa,
                                float outlierDeltaFactor, float dc, float kappa,
                                const unsigned int num_elements);
    DECLARE_TASKTYPE_AND_KERNEL(TAcc, AssignClusters,
                                unsigned int *numberOfClustersScalar,
                                bool writeOutNumClusters = true);
    DeviceRawPointers ptrs_;
  };

  CLUEAlgoAlpaka(TQueue &queue, float dc, float kappa, float outlierDeltaFactor,
                 bool verbose)
      : CLUEAlgo<T, NLAYERS>(dc, kappa, outlierDeltaFactor, verbose),
        device_(alpaka::getDev(queue)), queue_(queue),
        host_(
            alpaka::getDevByIdx(alpaka::Platform<alpaka::PlatformCpu>{}, 0u)) {}

  CLUEAlgoAlpaka(float dc, float kappa, float outlierDeltaFactor, bool verbose,
                 bool useAbsoluteSigma = false)
      : CLUEAlgo<T, NLAYERS>(dc, kappa, outlierDeltaFactor, verbose,
                             useAbsoluteSigma),
        device_(alpaka::getDevByIdx(alpaka::Platform<TAcc>{}, 0u)),
        queue_(device_), host_(alpaka::getDevByIdx(
                             alpaka::Platform<alpaka::PlatformCpu>{}, 0u)) {
    init_device();
  }

  ~CLUEAlgoAlpaka() { free_device(); }

  void makeClusters();
  void makeClustersCMSSW(const unsigned int points, const float *x,
                         const float *y, const int *layer, const float *weight,
                         const float *sigmaNoise, const uint32_t *detid,
                         float *rho, float *delta, unsigned int *nearestHigher,
                         int *clusterIndex, uint8_t *isSeed,
                         unsigned int *numberOfClustersScalar);

  // Device runner to submit kernels
  DeviceRunner device_runner_;

private:
  alpaka::Dev<TAcc> device_;
  // choose between Blocking and NonBlocking
  TQueue queue_;
  alpaka::DevCpu host_;

  // Memory management variables
  PointsBuf device_bufs_;
  std::optional<BufLayerTiles> device_hist_;
  std::optional<BufVecArrSeeds> device_seeds_;
  std::optional<BufVecArrFollowers> device_followers_;

  void init_device() {
    Idx const reserve = 1000000;
    // If Dim is not 1, fail compilation. This is assumed to be a
    // mono-dimensional problem
    static_assert(Dim::value == 1u);
    alpaka::Vec<Dim, Idx> const extents(reserve);

    // INPUT VARIABLES
    // Allocate device memory
    device_bufs_.x =
        std::make_optional(alpaka::allocBuf<float, Idx>(device_, extents));
    device_bufs_.y =
        std::make_optional(alpaka::allocBuf<float, Idx>(device_, extents));
    device_bufs_.layer =
        std::make_optional(alpaka::allocBuf<int, Idx>(device_, extents));
    device_bufs_.weight =
        std::make_optional(alpaka::allocBuf<float, Idx>(device_, extents));
    device_bufs_.sigmaNoise =
        std::make_optional(alpaka::allocBuf<float, Idx>(device_, extents));

    // RESULT VARIABLES
    device_bufs_.rho =
        std::make_optional(alpaka::allocBuf<float, Idx>(device_, extents));
    device_bufs_.delta =
        std::make_optional(alpaka::allocBuf<float, Idx>(device_, extents));
    device_bufs_.nearestHigher = std::make_optional(
        alpaka::allocBuf<unsigned int, Idx>(device_, extents));
    device_bufs_.clusterIndex =
        std::make_optional(alpaka::allocBuf<int, Idx>(device_, extents));
    device_bufs_.isSeed =
        std::make_optional(alpaka::allocBuf<uint8_t, Idx>(device_, extents));

    // INTERNAL VARIABLES
    alpaka::Vec<Dim, Idx> const layerTilesExtents(static_cast<Idx>(NLAYERS));
    device_hist_ = std::make_optional(
        alpaka::allocBuf<LayerTilesAcc, Idx>(device_, layerTilesExtents));

    alpaka::Vec<Dim, Idx> const seedsExtents(1u);
    device_seeds_ = std::make_optional(
        alpaka::allocBuf<GPUAlpaka::VecArray<int, maxNSeeds>, Idx>(
            device_, seedsExtents));

    device_followers_ = std::make_optional(
        alpaka::allocBuf<GPUAlpaka::VecArray<int, maxNFollowers>, Idx>(
            device_, extents));

    // Update RAW device pointers, grouped in a struct for convenience
    device_runner_.ptrs_.x = alpaka::getPtrNative(device_bufs_.x.value());
    device_runner_.ptrs_.y = alpaka::getPtrNative(device_bufs_.y.value());
    device_runner_.ptrs_.layer =
        alpaka::getPtrNative(device_bufs_.layer.value());
    device_runner_.ptrs_.weight =
        alpaka::getPtrNative(device_bufs_.weight.value());
    if (useAbsoluteSigma_)
      device_runner_.ptrs_.sigmaNoise =
          alpaka::getPtrNative(device_bufs_.sigmaNoise.value());

    // RESULT VARIABLES
    device_runner_.ptrs_.rho = alpaka::getPtrNative(device_bufs_.rho.value());
    device_runner_.ptrs_.delta =
        alpaka::getPtrNative(device_bufs_.delta.value());
    device_runner_.ptrs_.nearestHigher =
        alpaka::getPtrNative(device_bufs_.nearestHigher.value());
    device_runner_.ptrs_.clusterIndex =
        alpaka::getPtrNative(device_bufs_.clusterIndex.value());
    device_runner_.ptrs_.isSeed =
        alpaka::getPtrNative(device_bufs_.isSeed.value());

    // UPDATE RAW POINTERS FOR INTERNATL DATA STRUCTURES
    device_runner_.ptrs_.hist_ = alpaka::getPtrNative(device_hist_.value());
    device_runner_.ptrs_.seeds_ = alpaka::getPtrNative(device_seeds_.value());
    device_runner_.ptrs_.followers_ =
        alpaka::getPtrNative(device_followers_.value());
  }

  void free_device() {
    // Nothing really to be done here, since Alpaka memory buffers are
    // reference counted and will be automatically deleted when needed.
  }

  // Returns a view host using the memory allocated by the passed in
  // std::vector. The size of the view is inferred from the size of the vector.
  // This means the view will, possibly, become invalid if, in the meantime,
  // the vector re-allocated its underlying storage.
  template <typename TT>
  auto getViewHost(TT &t) -> ViewHostT<typename TT::value_type> {
    using type = typename TT::value_type;
    using Dim1 = alpaka::DimInt<1ul>;
    alpaka::Vec<Dim1, Idx> vectorSize(static_cast<Idx>(t.size()));
    ViewHostT<type> tempHostView(t.data(), host_, vectorSize);
    return tempHostView;
  }

  template <typename TT> auto getViewHost(TT *t, int size) -> ViewHostT<TT> {
    using Dim1 = alpaka::DimInt<1ul>;
    alpaka::Vec<Dim1, Idx> vectorSize(static_cast<Idx>(size));
    ViewHostT<TT> tempHostView(t, host_, vectorSize);
    return tempHostView;
  }

  void copy_todevice() {
    // input variables
    using Dim1 = alpaka::DimInt<1ul>;
    alpaka::Vec<Dim1, Idx> const extentToTransfer(static_cast<Idx>(points_.n));
    alpaka::memcpy(queue_, device_bufs_.x.value(),
                   getViewHost(points_.p_x, points_.n), extentToTransfer);
    alpaka::memcpy(queue_, device_bufs_.y.value(),
                   getViewHost(points_.p_y, points_.n), extentToTransfer);
    alpaka::memcpy(queue_, device_bufs_.layer.value(),
                   getViewHost(points_.p_layer, points_.n), extentToTransfer);
    alpaka::memcpy(queue_, device_bufs_.weight.value(),
                   getViewHost(points_.p_weight, points_.n), extentToTransfer);
    if (useAbsoluteSigma_)
      alpaka::memcpy(queue_, device_bufs_.sigmaNoise.value(),
                     getViewHost(points_.p_sigmaNoise, points_.n),
                     extentToTransfer);
    alpaka::wait(queue_);
  }

  void clear_internal_buffers() {
    // result variables
    using Dim1 = alpaka::DimInt<1ul>;
    alpaka::Vec<Dim1, Idx> extents(static_cast<Idx>(points_.n));
    alpaka::memset(queue_, device_bufs_.rho.value(), 0x0, extents);
    alpaka::memset(queue_, device_bufs_.delta.value(), 0x0, extents);
    alpaka::memset(queue_, device_bufs_.nearestHigher.value(), 0x0, extents);
    alpaka::memset(queue_, device_bufs_.clusterIndex.value(), 0x0, extents);
    alpaka::memset(queue_, device_bufs_.isSeed.value(), 0x0, extents);
    // algorithm internal variables
    // INTERNAL VARIABLES
    alpaka::Vec<Dim, Idx> const layerTilesExtents(static_cast<Idx>(NLAYERS));
    alpaka::memset(queue_, device_hist_.value(), 0x0, layerTilesExtents);

    alpaka::Vec<Dim, Idx> const seedsExtents(1u);
    alpaka::memset(queue_, device_seeds_.value(), 0x0, seedsExtents);

    alpaka::memset(queue_, device_followers_.value(), 0x0, extents);
    alpaka::wait(queue_);
  }

  void copy_tohost() {
    // result variables
    using Dim1 = alpaka::DimInt<1ul>;
    alpaka::Vec<Dim1, Idx> extents(static_cast<Idx>(points_.n));

    auto clusterHV = getViewHost(points_.clusterIndex);
    alpaka::memcpy(queue_, clusterHV, device_bufs_.clusterIndex.value(),
                   extents);
    if (verbose_) {
      // other variables, copy only when verbose_==True
      auto rhoHV = getViewHost(points_.rho);
      alpaka::memcpy(queue_, rhoHV, device_bufs_.rho.value(), extents);

      auto deltaHV = getViewHost(points_.delta);
      alpaka::memcpy(queue_, deltaHV, device_bufs_.delta.value(), extents);

      auto nearestHV = getViewHost(points_.nearestHigher);
      alpaka::memcpy(queue_, nearestHV, device_bufs_.nearestHigher.value(),
                     extents);

      auto isSeedHV = getViewHost(points_.isSeed);
      alpaka::memcpy(queue_, isSeedHV, device_bufs_.isSeed.value(), extents);
    }
    alpaka::wait(queue_);
  }
};

template <typename TAcc, typename TQueue, typename T, int NLAYERS>
ALPAKA_FN_ACC auto
CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::operator()(
    TAcc const &acc,
    CLUEAlgoAlpaka<TAcc, TQueue, T,
                   NLAYERS>::DeviceRunner::KernelComputeHistogram dummy,
    const unsigned int numberOfPoints) const -> void {
  const Idx i(alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u]);
  if (i < numberOfPoints) {
    // push index of points into tiles
    ptrs_.hist_[ptrs_.layer[i]].fill(ptrs_.x[i], ptrs_.y[i], i);
  }
}

template <typename TAcc, typename TQueue, typename T, int NLAYERS>
ALPAKA_FN_ACC auto
CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::operator()(
    TAcc const &acc,
    CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::KernelSortHistogram
        dummy) const -> void {
  /*
  const Idx layer(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
  ptrs_.hist_[layer].sort_unsafe();
  */
  const Idx i(alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u]);
  if (i < NLAYERS * T::nTiles) {
    int layer = i / T::nTiles;
    int bin = i - layer * T::nTiles;
    ptrs_.hist_[layer].sort_unsafe(bin);
  }
}

template <typename TAcc, typename TQueue, typename T, int NLAYERS>
ALPAKA_FN_ACC auto
CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::operator()(
    TAcc const &acc,
    CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::KernelInitializeHistogram
        dummy, const unsigned int numberOfPoints, const float outlierDeltaFactor, const float dc) const -> void {
  // Reset only the Tiles we will be effectively using.
  const Idx i(alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u]);
  if (i < numberOfPoints) {
    const float dm = outlierDeltaFactor * dc;
    float xi = ptrs_.x[i];
    float yi = ptrs_.y[i];
    const int layeri = ptrs_.layer[i];
    // get search box
    int4 search_box =
      ptrs_.hist_[layeri].searchBox(xi - dm, xi + dm, yi - dm, yi + dm);

    // loop over bins in the search box
    for (int xBin = search_box.x; xBin < search_box.y + 1; ++xBin) {
      for (int yBin = search_box.z; yBin < search_box.w + 1; ++yBin) {
        // get the id of this bin
        int globalBin = ptrs_.hist_[layeri].getGlobalBinByBin(xBin, yBin);
        ptrs_.hist_[layeri][globalBin].reset();
      }
    }
  }
}

template <typename TAcc, typename TQueue, typename T, int NLAYERS>
ALPAKA_FN_ACC auto
CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::operator()(
    TAcc const &acc,
    CLUEAlgoAlpaka<TAcc, TQueue, T,
                   NLAYERS>::DeviceRunner::KernelComputeLocalDensity dummy,
    float dc, const unsigned int numberOfPoints) const -> void {
  const Idx i(alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u]);
  printf("%s %i\n", __FILE__, __LINE__);
  if (i < numberOfPoints) {
    float rhoi{0.};
    int layeri = ptrs_.layer[i];
    float xi = ptrs_.x[i];
    float yi = ptrs_.y[i];

    // get search box
    int4 search_box =
        ptrs_.hist_[layeri].searchBox(xi - dc, xi + dc, yi - dc, yi + dc);

    // loop over bins in the search box
    for (int xBin = search_box.x; xBin < search_box.y + 1; ++xBin) {
      for (int yBin = search_box.z; yBin < search_box.w + 1; ++yBin) {
        // get the id of this bin
        int binId = ptrs_.hist_[layeri].getGlobalBinByBin(xBin, yBin);
        // get the size of this bin
        printf("BinId %i | (xBin, yBin) = (%i, %i) | ", binId, xBin, yBin);
        int binSize = ptrs_.hist_[layeri][binId].size();
        printf("BinSize %i\n", binSize);
        if (static_cast<unsigned int>(binSize) > numberOfPoints) continue;
        // interate inside this bin
        for (int binIter = 0; binIter < binSize; binIter++) {
          unsigned int j = ptrs_.hist_[layeri][binId][binIter];
          // query N_{dc_}(i)
          float xj = ptrs_.x[j];
          float yj = ptrs_.y[j];
          float dist_ij = ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj));
          if (dist_ij < dc * dc) {
            // sum weights within N_{dc_}(i)
            rhoi += (i == j ? 1.f : 0.5f) * ptrs_.weight[j];
          }
        } // end of interate inside this bin
      }
    } // end of loop over bins in search box
    ptrs_.rho[i] = rhoi;
  }
}

template <typename TAcc, typename TQueue, typename T, int NLAYERS>
ALPAKA_FN_ACC auto
CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::operator()(
    TAcc const &acc,
    CLUEAlgoAlpaka<TAcc, TQueue, T,
                   NLAYERS>::DeviceRunner::KernelComputeDistanceToHigherNoDetId
        dummy,
    float outlierDeltaFactor, float dc, const unsigned int numberOfPoints) const
    -> void {
  const Idx i(alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u]);
  float dm = outlierDeltaFactor * dc;

  if (i < numberOfPoints) {

    int layeri = ptrs_.layer[i];

    float deltai = std::numeric_limits<float>::max();
    unsigned int nearestHigheri = std::numeric_limits<unsigned int>::max();
    float xi = ptrs_.x[i];
    float yi = ptrs_.y[i];
    float rhoi = ptrs_.rho[i];
    float rho_max = 0.f;

    // get search box
    int4 search_box =
        ptrs_.hist_[layeri].searchBox(xi - dm, xi + dm, yi - dm, yi + dm);

    // loop over all bins in the search box
    for (int xBin = search_box.x; xBin < search_box.y + 1; ++xBin) {
      for (int yBin = search_box.z; yBin < search_box.w + 1; ++yBin) {
        // get the id of this bin
        int binId = ptrs_.hist_[layeri].getGlobalBinByBin(xBin, yBin);
        // get the size of this bin
        int binSize = ptrs_.hist_[layeri][binId].size();

        // iterate inside this bin
#if ORDER_TILE
        unsigned int old_j = 0;
#endif
        for (int binIter = 0; binIter < binSize; binIter++) {
          unsigned int j = ptrs_.hist_[layeri][binId][binIter];
#if ORDER_TILE
          assert(j >= old_j);
          old_j = j;
#endif
          // query N'_{dm}(i)
          float xj = ptrs_.x[j];
          float yj = ptrs_.y[j];
          float dist_ij =
              // std::sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj));
              ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj));
          bool is_inside = (dist_ij <= dm * dm);
          bool foundHigher = is_inside && (ptrs_.rho[j] > rhoi);
          // in the rare case where rho is the same, use indices
          foundHigher =
              is_inside && (foundHigher || ((ptrs_.rho[j] == rhoi) && (j > i)));
          if (foundHigher && dist_ij < deltai) {
            rho_max = ptrs_.rho[j];
            deltai = dist_ij;
            nearestHigheri = j;
          } else if (foundHigher && dist_ij == deltai &&
                     ptrs_.rho[j] > rho_max) {
            rho_max = ptrs_.rho[j];
            deltai = dist_ij;
            nearestHigheri = j;
          } else if (foundHigher && dist_ij == deltai &&
                     ptrs_.rho[j] == rho_max && j > i) {
            rho_max = ptrs_.rho[j];
            deltai = dist_ij;
            nearestHigheri = j;
          }
        } // end of interate inside this bin
      }
    } // end of loop over bins in search box
    ptrs_.delta[i] = std::sqrt(deltai);
    ptrs_.nearestHigher[i] = nearestHigheri;
  }
}

template <typename TAcc, typename TQueue, typename T, int NLAYERS>
ALPAKA_FN_ACC auto
CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::operator()(
    TAcc const &acc,
    CLUEAlgoAlpaka<TAcc, TQueue, T,
                   NLAYERS>::DeviceRunner::KernelComputeDistanceToHigher dummy,
    float outlierDeltaFactor, float dc, const unsigned int numberOfPoints) const
    -> void {
  const Idx i(alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u]);
  float dm = outlierDeltaFactor * dc;

  if (i < numberOfPoints) {

    int layeri = ptrs_.layer[i];

    float deltai = std::numeric_limits<float>::max();
    unsigned int nearestHigheri = std::numeric_limits<unsigned int>::max();
    float xi = ptrs_.x[i];
    float yi = ptrs_.y[i];
    float rhoi = ptrs_.rho[i];
    float rho_max = 0.f;

    // get search box
    int4 search_box =
        ptrs_.hist_[layeri].searchBox(xi - dm, xi + dm, yi - dm, yi + dm);

    // loop over all bins in the search box
    for (int xBin = search_box.x; xBin < search_box.y + 1; ++xBin) {
      for (int yBin = search_box.z; yBin < search_box.w + 1; ++yBin) {
        // get the id of this bin
        int binId = ptrs_.hist_[layeri].getGlobalBinByBin(xBin, yBin);
        // get the size of this bin
        int binSize = ptrs_.hist_[layeri][binId].size();

        // iterate inside this bin
#if ORDER_TILE
        unsigned int old_j = 0;
#endif
        for (int binIter = 0; binIter < binSize; binIter++) {
          unsigned int j = ptrs_.hist_[layeri][binId][binIter];
#if ORDER_TILE
          assert(j >= old_j);
          old_j = j;
#endif
          // query N'_{dm}(i)
          float xj = ptrs_.x[j];
          float yj = ptrs_.y[j];
          printf("X: %f, Y: %f\n", xj, yj);
          float dist_ij =
              // std::sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj));
              ((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj));
          bool foundHigher = (ptrs_.rho[j] > rhoi);
          // in the rare case where rho is the same, use detid
          foundHigher = foundHigher | ((ptrs_.rho[j] == rhoi) & (ptrs_.detid[j] > ptrs_.detid[i]));

          // combine the three conditions to avoid branching into different code segments to speedup the
          // memory write a little bit
          bool condition = foundHigher & ((dist_ij < deltai)
            | ( (dist_ij == deltai)
              & ((ptrs_.rho[j] > rho_max)
                | ((ptrs_.rho[j] == rho_max) & (ptrs_.detid[j] > ptrs_.detid[i])))));

          if(condition)
          {
            rho_max = ptrs_.rho[j];
            deltai = dist_ij;
            nearestHigheri = j;
          }
        } // end of interate inside this bin
      }
    } // end of loop over bins in search box
    ptrs_.delta[i] = std::sqrt(deltai);
    ptrs_.nearestHigher[i] = nearestHigheri;
  }
}

template <typename TAcc, typename TQueue, typename T, int NLAYERS>
ALPAKA_FN_ACC auto
CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::operator()(
    TAcc const &acc, KernelFindClustersKappa dummy, float outlierDeltaFactor,
    float dc, float kappa, const unsigned int numberOfPoints) const -> void {
  const Idx i(alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u]);

  if (i < numberOfPoints) {
    // initialize clusterIndex
    ptrs_.clusterIndex[i] = -1;
    // determine seed or outlier
    float deltai = ptrs_.delta[i];
    float rhoi = ptrs_.rho[i];
    float rhoc = ptrs_.sigmaNoise[i] * kappa;
    bool isSeed = (deltai > dc) && (rhoi >= rhoc);
    bool isOutlier = (deltai > outlierDeltaFactor * dc) && (rhoi < rhoc);

    if (isSeed) {
      // set isSeed as 1
      ptrs_.isSeed[i] = 1;
      ptrs_.seeds_[0].push_back(acc, i); // head of device_seeds_
    } else {
      ptrs_.isSeed[i] = 0;
      if (!isOutlier) {
        assert(ptrs_.nearestHigher[i] < numberOfPoints);
        // register as follower at its nearest higher
        ptrs_.followers_[ptrs_.nearestHigher[i]].push_back(acc, i);
      }
    }
  }
}

template <typename TAcc, typename TQueue, typename T, int NLAYERS>
ALPAKA_FN_ACC auto
CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::operator()(
    TAcc const &acc, KernelFindClusters dummy, float outlierDeltaFactor,
    float dc, float rhoc, const unsigned int numberOfPoints) const -> void {
  const Idx i(alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u]);

  if (i < numberOfPoints) {
    // initialize clusterIndex
    ptrs_.clusterIndex[i] = -1;
    // determine seed or outlier
    float deltai = ptrs_.delta[i];
    float rhoi = ptrs_.rho[i];
    bool isSeed = (deltai > dc) && (rhoi >= rhoc);
    bool isOutlier = (deltai > outlierDeltaFactor * dc) && (rhoi < rhoc);

    if (isSeed) {
      // set isSeed as 1
      ptrs_.isSeed[i] = 1;
      ptrs_.seeds_[0].push_back(acc, i); // head of device_seeds_
    } else {
      ptrs_.isSeed[i] = 0;
      if (!isOutlier) {
        assert(ptrs_.nearestHigher[i] < numberOfPoints);
        // register as follower at its nearest higher
        ptrs_.followers_[ptrs_.nearestHigher[i]].push_back(acc, i);
      }
    }
  }
}

template <typename TAcc, typename TQueue, typename T, int NLAYERS>
ALPAKA_FN_ACC auto
CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::operator()(
    TAcc const &acc,
    CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::KernelAssignClusters
        dummy,
    unsigned int *numberOfClustersScalar, bool writeOutNumClusters) const
    -> void {
  const Idx idxCls(alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u]);
  if (idxCls == 0 && writeOutNumClusters) {
    *numberOfClustersScalar = ptrs_.seeds_[0].size();
  }
  if (idxCls < (unsigned int)ptrs_.seeds_[0].size()) {
    int localStack[localStackSizePerSeed] = {-1};
    int localStackSize = 0;

    // assign cluster to seed[idxCls]
    int idxThisSeed = ptrs_.seeds_[0][idxCls];
    ptrs_.clusterIndex[idxThisSeed] = idxCls;
    // push_back idThisSeed to localStack
    assert(localStackSize < localStackSizePerSeed);
    localStack[localStackSize] = idxThisSeed;
    localStackSize++;

    // process all elements in localStack
    while (localStackSize > 0) {
      // get last element of localStack
      assert(localStackSize - 1 < localStackSizePerSeed);
      int idxEndOflocalStack = localStack[localStackSize - 1];

      int temp_clusterIndex = ptrs_.clusterIndex[idxEndOflocalStack];
      // pop_back last element of localStack
      assert(localStackSize - 1 < localStackSizePerSeed);
      localStack[localStackSize - 1] = -1;
      localStackSize--;

      // loop over followers of last element of localStack
      for (int j : ptrs_.followers_[idxEndOflocalStack]) {
        // pass id to follower
        ptrs_.clusterIndex[j] = temp_clusterIndex;
        // push_back follower to localStack
        assert(localStackSize < localStackSizePerSeed);
        localStack[localStackSize] = j;
        localStackSize++;
      }
    }
  }
}

template <typename TAcc, typename TQueue, typename T, int NLAYERS>
void CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::makeClusters() {
  copy_todevice();
  clear_internal_buffers();

  // Dimension the grid for submission
  alpaka::Vec<Dim, Idx> const threadsPerBlock(1024u);
  alpaka::Vec<Dim, Idx> const blocksPerGrid(
      static_cast<Idx>(ceil(points_.n / (float)threadsPerBlock[0])));
  alpaka::Vec<Dim, Idx> const elementsPerThread(1u);
  using WorkDiv = alpaka::WorkDivMembers<Dim, Idx>;
  auto const manualWorkDiv =
      WorkDiv{blocksPerGrid, threadsPerBlock, elementsPerThread};

  // Create the kernel execution tasks.
  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelComputeHistogram
      taskComputeHistogram;
  auto const kernelComputeHistogram = (alpaka::createTaskKernel<TAcc>(
      manualWorkDiv, device_runner_, taskComputeHistogram,
      static_cast<int>(points_.n)));

  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelComputeLocalDensity
      taskComputeLocalDensity;
  auto const kernelComputeLocalDensity = (alpaka::createTaskKernel<TAcc>(
      manualWorkDiv, device_runner_, taskComputeLocalDensity, dc_,
      static_cast<int>(points_.n)));

  typename CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::DeviceRunner::
      KernelComputeDistanceToHigherNoDetId taskComputeDistanceToHigherNoDetId;
  auto const kernelComputeDistanceToHigherNoDetId =
      (alpaka::createTaskKernel<TAcc>(
          manualWorkDiv, device_runner_, taskComputeDistanceToHigherNoDetId,
          outlierDeltaFactor_, dc_, static_cast<int>(points_.n)));

  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelFindClusters
      taskFindClusters;
  auto const kernelFindClusters = (alpaka::createTaskKernel<TAcc>(
      manualWorkDiv, device_runner_, taskFindClusters, outlierDeltaFactor_, dc_,
      rhoc_, static_cast<int>(points_.n)));

  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelFindClustersKappa
      taskFindClustersKappa;
  auto const kernelFindClustersKappa = (alpaka::createTaskKernel<TAcc>(
      manualWorkDiv, device_runner_, taskFindClustersKappa, outlierDeltaFactor_,
      dc_, kappa_, static_cast<int>(points_.n)));

  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelAssignClusters
      taskAssignClusters;
  auto const kernelAssignClusters = (alpaka::createTaskKernel<TAcc>(
      manualWorkDiv, device_runner_, taskAssignClusters, nullptr, false));

  // Enqueue the kernel execution task
  auto start = std::chrono::high_resolution_clock::now();
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed;

  start = std::chrono::high_resolution_clock::now();
  alpaka::enqueue(queue_, kernelComputeHistogram);
  alpaka::wait(queue_); // wait in case we are using an asynchronous queue to
  // time actual kernel runtime
  finish = std::chrono::high_resolution_clock::now();
  elapsed = finish - start;
  std::cout << "--- computeHistogram:            " << elapsed.count() * 1000
            << "ms\n";

  start = std::chrono::high_resolution_clock::now();
  alpaka::enqueue(queue_, kernelComputeLocalDensity);
  alpaka::wait(queue_); // wait in case we are using an asynchronous queue to
  // time actual kernel runtime
  finish = std::chrono::high_resolution_clock::now();
  elapsed = finish - start;
  std::cout << "--- computeLocalDensity:            " << elapsed.count() * 1000
            << "ms\n";

  start = std::chrono::high_resolution_clock::now();
  alpaka::enqueue(queue_, kernelComputeDistanceToHigherNoDetId);
  alpaka::wait(queue_); // wait in case we are using an asynchronous queue to
  // time actual kernel runtime
  finish = std::chrono::high_resolution_clock::now();
  elapsed = finish - start;
  std::cout << "--- computeDistanceToHigher:            "
            << elapsed.count() * 1000 << "ms\n";

  start = std::chrono::high_resolution_clock::now();
  if (useAbsoluteSigma_) {
    alpaka::enqueue(queue_, kernelFindClustersKappa);
  } else {
    alpaka::enqueue(queue_, kernelFindClusters);
  }
  alpaka::wait(queue_); // wait in case we are using an asynchronous queue to
  // time actual kernel runtime
  finish = std::chrono::high_resolution_clock::now();
  elapsed = finish - start;
  std::cout << "--- findClusters:            " << elapsed.count() * 1000
            << "ms\n";

  start = std::chrono::high_resolution_clock::now();
  alpaka::enqueue(queue_, kernelAssignClusters);
  alpaka::wait(queue_); // wait in case we are using an asynchronous queue to
  // time actual kernel runtime
  finish = std::chrono::high_resolution_clock::now();
  elapsed = finish - start;
  std::cout << "--- assignClusters:            " << elapsed.count() * 1000
            << "ms\n";

  copy_tohost();
}

template <typename TAcc, typename TQueue, typename T, int NLAYERS>
void CLUEAlgoAlpaka<TAcc, TQueue, T, NLAYERS>::makeClustersCMSSW(
    const unsigned int points, const float *x, const float *y, const int *layer,
    const float *weight, const float *sigmaNoise, const uint32_t *detid,
    float *rho, float *delta, unsigned int *nearestHigher, int *clusterIndex,
    uint8_t *isSeed, unsigned int *numberOfClustersScalar) {

  // std::cout << "makeClustersCMSSW received " << points << " RecHits" <<
  // std::endl;
  using WorkDiv = alpaka::WorkDivMembers<Dim, Idx>;

  // INTERNAL TILES VARIABLES
  Idx const reserve = points;
  // If Dim is not 1, fail compilation. This is assumed to be a
  // mono-dimensional problem
  static_assert(Dim::value == 1u);
  alpaka::Vec<Dim, Idx> const extents(reserve);

  alpaka::Vec<Dim, Idx> const layerTilesExtents(static_cast<Idx>(NLAYERS));
  device_hist_ = std::make_optional(
      alpaka::allocAsyncBuf<LayerTilesAcc, Idx>(queue_, layerTilesExtents));
  alpaka::Vec<Dim, Idx> const seedsExtents(1u);
  device_seeds_ = std::make_optional(
      alpaka::allocAsyncBuf<GPUAlpaka::VecArray<int, maxNSeeds>, Idx>(
          queue_, seedsExtents));
  device_followers_ = std::make_optional(
      alpaka::allocAsyncBuf<GPUAlpaka::VecArray<int, maxNFollowers>, Idx>(
          queue_, extents));
  
  auto numberOfClustersView =
    alpaka::createView(alpaka::getDev(queue_), numberOfClustersScalar, 1u);
  alpaka::memset(queue_, numberOfClustersView, 0x0);
  alpaka::memset(queue_, device_seeds_.value(), 0x0, seedsExtents);
  alpaka::memset(queue_, device_followers_.value(), 0x0, extents);

  // Set Device Raw Pointers using values from outsice and also internal buffers
  device_runner_.ptrs_.x = const_cast<float *>(x);
  device_runner_.ptrs_.y = const_cast<float *>(y);
  device_runner_.ptrs_.layer = const_cast<int *>(layer);
  device_runner_.ptrs_.weight = const_cast<float *>(weight);
  device_runner_.ptrs_.sigmaNoise = const_cast<float *>(sigmaNoise);
  device_runner_.ptrs_.detid = const_cast<uint32_t *>(detid);

  // RESULT VARIABLES
  device_runner_.ptrs_.rho = rho;
  device_runner_.ptrs_.delta = delta;
  device_runner_.ptrs_.nearestHigher = nearestHigher;
  device_runner_.ptrs_.clusterIndex = clusterIndex;
  device_runner_.ptrs_.isSeed = isSeed;

  // UPDATE RAW POINTERS FOR INTERNATL DATA STRUCTURES
  device_runner_.ptrs_.hist_ = alpaka::getPtrNative(device_hist_.value());
  device_runner_.ptrs_.seeds_ = alpaka::getPtrNative(device_seeds_.value());
  device_runner_.ptrs_.followers_ =
      alpaka::getPtrNative(device_followers_.value());


  // Dimension the grid for submission
  Idx threads_per_block = 256u;
  if constexpr (std::is_same_v<alpaka::Dev<TAcc>, alpaka::DevCpu>) {
    threads_per_block = 1u;
  }
  alpaka::Vec<Dim, Idx> const threadsPerBlock(threads_per_block);

  alpaka::Vec<Dim, Idx> const blocksPerGrid(
      static_cast<Idx>(ceil(points / (float)threadsPerBlock[0])));
  alpaka::Vec<Dim, Idx> const elementsPerThread(1u);
  auto const manualWorkDiv =
      WorkDiv{blocksPerGrid, threadsPerBlock, elementsPerThread};

  // Set internal size of each bin in the tiles at 0, leave the rest unchanged
  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelInitializeHistogram
      taskInitializeHistogram;
  auto const kernelInitializeHistogram = alpaka::createTaskKernel<TAcc>(
      manualWorkDiv, device_runner_, taskInitializeHistogram, points, outlierDeltaFactor_, dc_);

  // Create the kernel execution tasks.
  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelComputeHistogram
      taskComputeHistogram;
  auto const kernelComputeHistogram = alpaka::createTaskKernel<TAcc>(
      manualWorkDiv, device_runner_, taskComputeHistogram, points);

#if ORDER_TILE
  // printf("Sorting tile for all layers.\n");
  alpaka::Vec<Dim, Idx> const threadsPerBlockTile(512u);
  alpaka::Vec<Dim, Idx> const blocksPerGridTile(
      std::ceil((NLAYERS * T::nTiles) / (float)threadsPerBlockTile[0]));
  auto const manualWorkDivTile =
      WorkDiv{blocksPerGridTile, threadsPerBlockTile, elementsPerThread};

  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelSortHistogram
      taskSortHistogram;
  auto const kernelSortHistogram = alpaka::createTaskKernel<TAcc>(
      manualWorkDivTile, device_runner_, taskSortHistogram);
#endif

  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelComputeLocalDensity
      taskComputeLocalDensity;
  auto const kernelComputeLocalDensity = (alpaka::createTaskKernel<TAcc>(
      manualWorkDiv, device_runner_, taskComputeLocalDensity, dc_,
      static_cast<int>(points)));

  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelComputeDistanceToHigher
      taskComputeDistanceToHigher;
  auto const kernelComputeDistanceToHigher = (alpaka::createTaskKernel<TAcc>(
      manualWorkDiv, device_runner_, taskComputeDistanceToHigher,
      outlierDeltaFactor_, dc_, static_cast<int>(points)));

  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelFindClustersKappa
      taskFindClustersKappa;
  auto const kernelFindClustersKappa = (alpaka::createTaskKernel<TAcc>(
      manualWorkDiv, device_runner_, taskFindClustersKappa, outlierDeltaFactor_,
      dc_, kappa_, static_cast<int>(points)));

  typename CLUEAlgoAlpaka<TAcc, TQueue, T,
                          NLAYERS>::DeviceRunner::KernelAssignClusters
      taskAssignClusters;
  auto const kernelAssignClusters = (alpaka::createTaskKernel<TAcc>(
      manualWorkDiv, device_runner_, taskAssignClusters,
      numberOfClustersScalar));

  // Enqueue the kernel execution task
  printf("Hello hello\n");
  alpaka::enqueue(queue_, kernelInitializeHistogram);
  alpaka::enqueue(queue_, kernelComputeHistogram);
#if ORDER_TILE
  alpaka::enqueue(queue_, kernelSortHistogram);
#endif
  alpaka::enqueue(queue_, kernelComputeLocalDensity);
  alpaka::enqueue(queue_, kernelComputeDistanceToHigher);
  alpaka::enqueue(queue_, kernelFindClustersKappa);
  alpaka::enqueue(queue_, kernelAssignClusters);
}
