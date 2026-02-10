#ifndef TilesConstants_h
#define TilesConstants_h

namespace util {
static constexpr int32_t ceil(float num) {
  return (static_cast<float>(static_cast<int32_t>(num)) == num)
             ? static_cast<int32_t>(num)
             : static_cast<int32_t>(num) + ((num > 0) ? 1 : 0);
}
}  // namespace util

struct TilesConstants {
  static constexpr float tileSize = 10.f;
  static constexpr float minDim1 = -250.f;
  static constexpr float maxDim1 = 250.f;
  static constexpr float minDim2 = -250.f;
  static constexpr float maxDim2 = 250.f;
  static constexpr int nColumns = util::ceil((maxDim1 - minDim1) / tileSize);
  static constexpr int nRows = util::ceil((maxDim2 - minDim2) / tileSize);
  static constexpr float invDim1BinSize = nColumns / (maxDim1 - minDim1);
  static constexpr float invDim2BinSize = nRows / (maxDim2 - minDim2);
  static constexpr int nTiles = nColumns * nRows;
  static constexpr int maxTileDepth = 64;  // For accelerators.
};

#endif  // TilesConstants_h
