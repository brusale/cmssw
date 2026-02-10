#ifndef GPUVecArrayAlpaka_h
#define GPUVecArrayAlpaka_h

namespace GPUAlpaka {

template <class T, int maxSize>
struct VecArray {
  inline constexpr int push_back_unsafe(const T &element) {
    auto previousSize = m_size;
    m_size++;
    if (previousSize < maxSize) {
      m_data[previousSize] = element;
      return previousSize;
    } else {
      --m_size;
      return -1;
    }
  }

  template <class... Ts>
  constexpr int emplace_back_unsafe(Ts &&...args) {
    auto previousSize = m_size;
    m_size++;
    if (previousSize < maxSize) {
      (new (&m_data[previousSize]) T(std::forward<Ts>(args)...));
      return previousSize;
    } else {
      --m_size;
      return -1;
    }
  }

  inline constexpr T &back() const {
    if (m_size > 0) {
      return m_data[m_size - 1];
    } else {
      assert(0);
      return T();  // undefined behaviour
    }
  }

  // thread-safe version of the vector, when used in a CUDA kernel
  template <typename T_Acc>
  ALPAKA_FN_ACC int push_back(const T_Acc &acc, const T &element) {
    auto previousSize = atomicAdd(acc, &m_size, 1, alpaka::hierarchy::Blocks{});
    if (previousSize < maxSize) {
      m_data[previousSize] = element;
      return previousSize;
    } else {
      assert(0);
      atomicSub(acc, &m_size, 1, alpaka::hierarchy::Blocks{});
      //assert(("Too few elemets reserved", maxSize));
      return -1;
    }
  }

  template <typename T_Acc, class... Ts>
  ALPAKA_FN_ACC int emplace_back(const T_Acc &acc, Ts &&...args) {
    auto previousSize = atomicAdd(acc, &m_size, 1, alpaka::hierarchy::Blocks{});
    if (previousSize < maxSize) {
      (new (&m_data[previousSize]) T(std::forward<Ts>(args)...));
      return previousSize;
    } else {
      assert(0);
      atomicSub(acc, &m_size, 1, alpaka::hierarchy::Blocks{});
      return -1;
    }
  }

  template <typename T_Acc, class... Ts>
  ALPAKA_FN_ACC inline T pop_back() {
    if (m_size > 0) {
      auto previousSize = m_size--;
      return m_data[previousSize - 1];
    } else {
      assert(0);
      return T();
    }
  }

  template <typename T_Acc>
    ALPAKA_FN_ACC inline void sort_unsafe(const T_Acc &acc) {
      for (int i = 0; i < m_size - 1; ++i) {
        for (int j = 0; j < m_size - i - 1; ++j) {
          if (m_data[j] > m_data[j + 1]) {
            // Swap elements if they are in the wrong order
            T temp = m_data[j];
            m_data[j] = m_data[j + 1];
            m_data[j + 1] = temp;
          }
        }
      }
    }

  inline constexpr T const *begin() const { return m_data; }
  inline constexpr T const *end() const { return m_data + m_size; }
  inline constexpr T *begin() { return m_data; }
  inline constexpr T *end() { return m_data + m_size; }
  inline constexpr int size() const { return m_size; }
  inline constexpr T &operator[](int i) { return m_data[i]; }
  inline constexpr const T &operator[](int i) const { return m_data[i]; }
  inline constexpr void reset() { m_size = 0; }
  inline constexpr int capacity() const { return maxSize; }
  inline constexpr T const *data() const { return m_data; }
  inline constexpr void resize(int size) { m_size = size; }
  inline constexpr bool empty() const { return 0 == m_size; }
  inline constexpr bool full() const { return maxSize == m_size; }

  int m_size = 0;

  T m_data[maxSize];
};

}  // namespace GPUAlpaka

#endif  // GPUVecArray_h
