#ifndef Points_h
#define Points_h

struct Points {
  // For internally managed input
  std::vector<float> x;
  std::vector<float> y;
  std::vector<int> layer;
  std::vector<float> weight;
  std::vector<float> sigmaNoise;

  // For externally managed input
  const float *p_x;
  const float *p_y;
  const int *p_layer;
  const float *p_weight;
  const float *p_sigmaNoise;

  std::vector<float> rho;
  std::vector<float> delta;
  std::vector<unsigned int> nearestHigher;
  std::vector<int> clusterIndex;
  std::vector<std::vector<int>> followers;
  std::vector<uint8_t> isSeed;

  int n;

  void clear() {
    x.clear();
    y.clear();
    layer.clear();
    weight.clear();
    sigmaNoise.clear();

    p_x = nullptr;
    p_y = nullptr;
    p_layer = nullptr;
    p_weight = nullptr;
    p_sigmaNoise = nullptr;

    rho.clear();
    delta.clear();
    nearestHigher.clear();
    clusterIndex.clear();
    followers.clear();
    isSeed.clear();

    n = 0;
  }
};
#endif
