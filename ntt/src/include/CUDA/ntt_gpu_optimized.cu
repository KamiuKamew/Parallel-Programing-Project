#include "ntt.h"
#include <cmath>
#include <iostream>
#include <vector>

#define CHECK_CUDA_OPTIMIZED(call)                                             \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      std::cerr << "CUDA error in " << __FILE__ << " at line " << __LINE__     \
                << ": " << cudaGetErrorString(err) << std::endl;               \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

// 引入基础模乘算法定义
enum MulMode { MUL_NAIVE = 0, MUL_MONT = 1, MUL_BARRETT = 2 };

// =============================================================================
// 基础GPU设备函数（从原实现复制）
// =============================================================================

template <typename T>
__device__ T gpu_mont_reduce_correct(__uint128_t t, T mod, T neg_r_inv) {
  constexpr int word_bits = sizeof(T) * 8;
  T m = (T)t * neg_r_inv;
  __uint128_t tmp = t + (__uint128_t)m * mod;
  T res = (T)(tmp >> word_bits);
  res = res - (mod & -(res >= mod));
  return res;
}

template <typename T>
__device__ T gpu_mont_add_correct(T a_mont, T b_mont, T mod) {
  return (a_mont + b_mont >= mod) ? (a_mont + b_mont - mod) : (a_mont + b_mont);
}

template <typename T>
__device__ T gpu_mont_sub_correct(T a_mont, T b_mont, T mod) {
  return (a_mont >= b_mont) ? (a_mont - b_mont) : (a_mont + mod - b_mont);
}

template <typename T>
__device__ T gpu_mont_mul_correct(T a_mont, T b_mont, T mod, T neg_r_inv) {
  return gpu_mont_reduce_correct((__uint128_t)a_mont * b_mont, mod, neg_r_inv);
}

__device__ inline u32 mod_mul_barrett_32(u32 a_mont, u32 b_mont, u32 p,
                                         u64 mu) {
  u32 inv = 1;
  for (int i = 0; i < 6; ++i) {
    inv = (u32)((__uint128_t)inv * (2 - (__uint128_t)p * inv));
  }
  u32 neg_r_inv = -inv;
  return gpu_mont_mul_correct(a_mont, b_mont, p, neg_r_inv);
}

template <int MODE, typename T>
__device__ inline T mod_mul_unified(T a, T b, T p, T neg_r_inv, u64 mu) {
  if constexpr (MODE == MUL_NAIVE) {
    return gpu_mont_mul_correct(a, b, p, neg_r_inv);
  } else if constexpr (MODE == MUL_MONT) {
    return gpu_mont_mul_correct(a, b, p, neg_r_inv);
  } else { // MUL_BARRETT
    if constexpr (sizeof(T) == 4) {
      return mod_mul_barrett_32(a, b, p, mu);
    } else {
      return gpu_mont_mul_correct(a, b, p, neg_r_inv);
    }
  }
}

// =============================================================================
// 优化策略1：共享内存优化的NTT核函数
// =============================================================================

template <typename T, int MODE>
__global__ void
ntt_stage_kernel_shared_memory_optimized(T *a, const T *twiddles, int mid,
                                         int n, T p, T neg_r_inv, u64 mu) {
  // 共享内存声明：预加载旋转因子
  extern __shared__ __align__(8) unsigned char shared_memory[];
  T *shared_twiddles = reinterpret_cast<T *>(shared_memory);

  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int local_tid = threadIdx.x;
  int total = n >> 1;

  if (tid >= total)
    return;

  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;
  int tw_idx = (mid - 1) + k;

  // 协作加载旋转因子到共享内存
  if (local_tid < mid && tw_idx < n) {
    shared_twiddles[local_tid] = twiddles[tw_idx];
  }
  __syncthreads();

  // 数据访问
  int data_idx1 = j + k;
  int data_idx2 = j + k + mid;

  T x = a[data_idx1];
  T y_data = a[data_idx2];

  // 使用共享内存中的旋转因子
  T w_mont =
      (k < blockDim.x && k < mid) ? shared_twiddles[k] : twiddles[tw_idx];

  T y = mod_mul_unified<MODE>(w_mont, y_data, p, neg_r_inv, mu);
  a[data_idx1] = gpu_mont_add_correct(x, y, p);
  a[data_idx2] = gpu_mont_sub_correct(x, y, p);
}

// 逆变换核函数（共享内存优化版本）
template <typename T, int MODE>
__global__ void ntt_stage_kernel_inv_shared_memory(T *a, const T *twiddles,
                                                   int mid, int n, T p,
                                                   T neg_r_inv, u64 mu) {
  extern __shared__ __align__(8) unsigned char shared_memory[];
  T *shared_twiddles = reinterpret_cast<T *>(shared_memory);

  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int local_tid = threadIdx.x;
  int total = n >> 1;

  if (tid >= total)
    return;

  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;
  int tw_idx = (mid - 1) + k;

  // 预加载旋转因子
  if (local_tid < mid && tw_idx < n) {
    shared_twiddles[local_tid] = twiddles[tw_idx];
  }
  __syncthreads();

  int data_idx1 = j + k;
  int data_idx2 = j + k + mid;

  T x = a[data_idx1];
  T y = a[data_idx2];

  a[data_idx1] = gpu_mont_add_correct(x, y, p);
  T temp = gpu_mont_sub_correct(x, y, p);

  T w_mont =
      (k < blockDim.x && k < mid) ? shared_twiddles[k] : twiddles[tw_idx];
  a[data_idx2] = mod_mul_unified<MODE>(w_mont, temp, p, neg_r_inv, mu);
}

// =============================================================================
// 优化策略2：向量化处理 - 每个线程处理多个元素
// =============================================================================

template <typename T, int MODE, int ELEMENTS_PER_THREAD>
__global__ void ntt_stage_kernel_vectorized(T *a, const T *twiddles, int mid,
                                            int n, T p, T neg_r_inv, u64 mu) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total = n >> 1;

// 每个线程处理多个蝶形操作
#pragma unroll
  for (int elem = 0; elem < ELEMENTS_PER_THREAD; elem++) {
    int current_tid = tid * ELEMENTS_PER_THREAD + elem;
    if (current_tid >= total)
      break;

    int step = mid << 1;
    int j = (current_tid / mid) * step;
    int k = current_tid % mid;
    int tw_idx = (mid - 1) + k;

    T w_mont = twiddles[tw_idx];
    T x = a[j + k];
    T y = mod_mul_unified<MODE>(w_mont, a[j + k + mid], p, neg_r_inv, mu);

    a[j + k] = gpu_mont_add_correct(x, y, p);
    a[j + k + mid] = gpu_mont_sub_correct(x, y, p);
  }
}

// =============================================================================
// 优化策略3：自适应线程配置选择
// =============================================================================

struct OptimalConfig {
  int block_size;
  int elements_per_thread;
  bool use_shared_memory;
  size_t shared_mem_size;
};

OptimalConfig select_optimal_config(int n, int gpu_sm_count) {
  OptimalConfig config;

  if (n <= 1024) {
    // 小规模：减少启动开销
    config.block_size = 128;
    config.elements_per_thread = 1;
    config.use_shared_memory = false;
    config.shared_mem_size = 0;
  } else if (n <= 16384) {
    // 中等规模：平衡内存和计算
    config.block_size = 256;
    config.elements_per_thread = 2;
    config.use_shared_memory = true;
    config.shared_mem_size = 256 * sizeof(u64);
  } else {
    // 大规模：最大化并行度
    config.block_size = 512;
    config.elements_per_thread = 4;
    config.use_shared_memory = true;
    config.shared_mem_size = 512 * sizeof(u64);
  }

  return config;
}

// =============================================================================
// 优化策略4：内存池管理
// =============================================================================

template <typename T> class GpuMemoryPool {
private:
  T *gpu_buffer_a;
  T *gpu_buffer_b;
  T *gpu_buffer_ab;
  T *gpu_twiddles;
  T *gpu_twiddles_inv;
  size_t current_capacity;

public:
  GpuMemoryPool() : current_capacity(0) {
    gpu_buffer_a = nullptr;
    gpu_buffer_b = nullptr;
    gpu_buffer_ab = nullptr;
    gpu_twiddles = nullptr;
    gpu_twiddles_inv = nullptr;
  }

  ~GpuMemoryPool() { cleanup(); }

  void ensure_capacity(size_t n_expanded) {
    if (n_expanded > current_capacity) {
      cleanup();

      CHECK_CUDA_OPTIMIZED(cudaMalloc(&gpu_buffer_a, n_expanded * sizeof(T)));
      CHECK_CUDA_OPTIMIZED(cudaMalloc(&gpu_buffer_b, n_expanded * sizeof(T)));
      CHECK_CUDA_OPTIMIZED(cudaMalloc(&gpu_buffer_ab, n_expanded * sizeof(T)));
      CHECK_CUDA_OPTIMIZED(cudaMalloc(&gpu_twiddles, n_expanded * sizeof(T)));
      CHECK_CUDA_OPTIMIZED(
          cudaMalloc(&gpu_twiddles_inv, n_expanded * sizeof(T)));

      current_capacity = n_expanded;
    }
  }

  void cleanup() {
    if (current_capacity > 0) {
      cudaFree(gpu_buffer_a);
      cudaFree(gpu_buffer_b);
      cudaFree(gpu_buffer_ab);
      cudaFree(gpu_twiddles);
      cudaFree(gpu_twiddles_inv);
      current_capacity = 0;
    }
  }

  T *get_buffer_a() { return gpu_buffer_a; }
  T *get_buffer_b() { return gpu_buffer_b; }
  T *get_buffer_ab() { return gpu_buffer_ab; }
  T *get_twiddles() { return gpu_twiddles; }
  T *get_twiddles_inv() { return gpu_twiddles_inv; }
};

// 全局内存池实例
static GpuMemoryPool<u32> gpu_pool_u32;
static GpuMemoryPool<u64> gpu_pool_u64;

// =============================================================================
// 基础算法函数（从原实现引用）
// =============================================================================

// 扩展数组大小到2的幂
template <typename T> T expand_n(T n) {
  T expanded = 1;
  while (expanded < n) {
    expanded <<= 1;
  }
  return expanded;
}

// 扩展数组内容
template <typename T> T *expand_a(T *a, T n, T n_expanded) {
  T *expanded = new T[n_expanded];
  for (T i = 0; i < n; i++) {
    expanded[i] = a[i];
  }
  for (T i = n; i < n_expanded; i++) {
    expanded[i] = 0;
  }
  return expanded;
}

// 位反转置换
template <typename T> void bit_reverse_permute(T *a, T n) {
  for (T i = 0; i < n; i++) {
    T j = 0;
    T temp = i;
    T log_n = 0;
    T temp_n = n;
    while (temp_n > 1) {
      temp_n >>= 1;
      log_n++;
    }

    for (T k = 0; k < log_n; k++) {
      j <<= 1;
      j |= (temp & 1);
      temp >>= 1;
    }

    if (i < j) {
      T swap = a[i];
      a[i] = a[j];
      a[j] = swap;
    }
  }
}

template <typename T>
void generate_twiddle_table(std::vector<T> &table, T n, T p, T omega,
                            bool inverse) {
  table.resize(n);
  MontMod<T> mont(p);
  T omega_mont = mont.from_T(omega);
  if (inverse)
    omega_mont = mont.inv(omega_mont);
  T one_mont = mont.from_T(1);

  for (int level = 0; (T(1) << level) < n; ++level) {
    T mid = T(1) << level;
    T exp = (p - 1) / (mid << 1);
    T Wn_mont = mont.pow(omega_mont, exp);
    T w = one_mont;
    T offset = mid - 1;
    for (T k = 0; k < mid; ++k) {
      table[offset + k] = w;
      w = mont.mul(w, Wn_mont);
    }
  }
}

template <typename T>
__global__ void pointwise_multiply_kernel_correct(T *a_mont, T *b_mont,
                                                  T *ab_mont, T n, T mod,
                                                  T neg_r_inv) {
  T idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < n) {
    ab_mont[idx] =
        gpu_mont_mul_correct(a_mont[idx], b_mont[idx], mod, neg_r_inv);
  }
}

// =============================================================================
// 优化的核心GPU NTT实现
// =============================================================================

template <typename T, int MODE>
void poly_multiply_ntt_gpu_optimized(T *a, T *b, T *ab, T n, T p, T omega) {
  MontMod<T> montMod(p);
  T n_expanded = expand_n(2 * n - 1);

  // 获取GPU硬件信息
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);

  OptimalConfig config =
      select_optimal_config(n_expanded, prop.multiProcessorCount);

  // 计算Montgomery参数
  T neg_r_inv;
  {
    constexpr int word_bits = sizeof(T) * 8;
    T inv = 1;
    int inv_iterations = (word_bits <= 64) ? 6 : 7;
    for (int i = 0; i < inv_iterations; ++i) {
      inv = (T)((__uint128_t)inv * (2 - (__uint128_t)p * inv));
    }
    neg_r_inv = -inv;
  }

  u64 mu = 0;
  if constexpr (MODE == MUL_BARRETT && sizeof(T) == 4) {
    mu = ((1ULL << 32) * (1ULL << 32)) / p;
  }

  // 使用内存池优化内存管理
  auto &gpu_pool = (sizeof(T) == 4)
                       ? reinterpret_cast<GpuMemoryPool<T> &>(gpu_pool_u32)
                       : reinterpret_cast<GpuMemoryPool<T> &>(gpu_pool_u64);

  gpu_pool.ensure_capacity(n_expanded);

  T *a_mont_gpu = gpu_pool.get_buffer_a();
  T *b_mont_gpu = gpu_pool.get_buffer_b();
  T *ab_mont_gpu = gpu_pool.get_buffer_ab();
  T *tw_gpu = gpu_pool.get_twiddles();
  T *tw_inv_gpu = gpu_pool.get_twiddles_inv();

  // 数据预处理
  T *a_expanded = expand_a(a, n, n_expanded);
  T *b_expanded = expand_a(b, n, n_expanded);
  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  // 预计算旋转因子
  std::vector<T> tw_host(n_expanded);
  std::vector<T> tw_inv_host(n_expanded);
  generate_twiddle_table(tw_host, n_expanded, p, omega, false);
  generate_twiddle_table(tw_inv_host, n_expanded, p, omega, true);

  // Host端Montgomery转换
  T *a_mont_host = new T[n_expanded];
  T *b_mont_host = new T[n_expanded];

  for (T i = 0; i < n_expanded; i++) {
    a_mont_host[i] = montMod.from_T(a_expanded[i]);
    b_mont_host[i] = montMod.from_T(b_expanded[i]);
  }

  // 异步数据传输
  cudaStream_t stream1, stream2;
  cudaStreamCreate(&stream1);
  cudaStreamCreate(&stream2);

  CHECK_CUDA_OPTIMIZED(cudaMemcpyAsync(a_mont_gpu, a_mont_host,
                                       n_expanded * sizeof(T),
                                       cudaMemcpyHostToDevice, stream1));
  CHECK_CUDA_OPTIMIZED(cudaMemcpyAsync(b_mont_gpu, b_mont_host,
                                       n_expanded * sizeof(T),
                                       cudaMemcpyHostToDevice, stream2));
  CHECK_CUDA_OPTIMIZED(cudaMemcpyAsync(tw_gpu, tw_host.data(),
                                       n_expanded * sizeof(T),
                                       cudaMemcpyHostToDevice, stream1));
  CHECK_CUDA_OPTIMIZED(cudaMemcpyAsync(tw_inv_gpu, tw_inv_host.data(),
                                       n_expanded * sizeof(T),
                                       cudaMemcpyHostToDevice, stream2));

  cudaStreamSynchronize(stream1);
  cudaStreamSynchronize(stream2);

  // 优化的NTT前向变换
  for (int level = 0; (T(1) << level) < n_expanded; ++level) {
    T mid = T(1) << level;
    int total_butterfly = n_expanded >> 1;

    // 根据配置选择最优核函数
    if (config.use_shared_memory && mid <= config.block_size) {
      int GRID = (total_butterfly + config.block_size - 1) / config.block_size;
      size_t shared_mem_size =
          std::min((size_t)mid, (size_t)config.block_size) * sizeof(T);

      ntt_stage_kernel_shared_memory_optimized<T, MODE>
          <<<GRID, config.block_size, shared_mem_size>>>(
              a_mont_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
      ntt_stage_kernel_shared_memory_optimized<T, MODE>
          <<<GRID, config.block_size, shared_mem_size>>>(
              b_mont_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
    } else if (config.elements_per_thread > 1) {
      int effective_threads =
          (total_butterfly + config.elements_per_thread - 1) /
          config.elements_per_thread;
      int GRID =
          (effective_threads + config.block_size - 1) / config.block_size;

      if (config.elements_per_thread == 2) {
        ntt_stage_kernel_vectorized<T, MODE, 2><<<GRID, config.block_size>>>(
            a_mont_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
        ntt_stage_kernel_vectorized<T, MODE, 2><<<GRID, config.block_size>>>(
            b_mont_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
      } else if (config.elements_per_thread == 4) {
        ntt_stage_kernel_vectorized<T, MODE, 4><<<GRID, config.block_size>>>(
            a_mont_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
        ntt_stage_kernel_vectorized<T, MODE, 4><<<GRID, config.block_size>>>(
            b_mont_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
      }
    }
    CHECK_CUDA_OPTIMIZED(cudaGetLastError());
  }
  CHECK_CUDA_OPTIMIZED(cudaDeviceSynchronize());

  // 点乘
  {
    dim3 blockSize(config.block_size);
    dim3 gridSize((n_expanded + blockSize.x - 1) / blockSize.x);
    pointwise_multiply_kernel_correct<<<gridSize, blockSize>>>(
        a_mont_gpu, b_mont_gpu, ab_mont_gpu, n_expanded, p, neg_r_inv);
    CHECK_CUDA_OPTIMIZED(cudaGetLastError());
  }
  CHECK_CUDA_OPTIMIZED(cudaDeviceSynchronize());

  // NTT逆变换（使用基础版本保证正确性）
  for (int level = (int)std::log2(n_expanded) - 1; level >= 0; --level) {
    T mid = T(1) << level;
    int total_butterfly = n_expanded >> 1;
    int GRID = (total_butterfly + config.block_size - 1) / config.block_size;

    // 使用基础逆变换核函数（保证正确性）
    ntt_stage_kernel_inv_shared_memory<T, MODE><<<GRID, config.block_size>>>(
        ab_mont_gpu, tw_inv_gpu, mid, n_expanded, p, neg_r_inv, mu);
    CHECK_CUDA_OPTIMIZED(cudaGetLastError());
  }
  CHECK_CUDA_OPTIMIZED(cudaDeviceSynchronize());

  // 结果处理
  T *temp_result = new T[n_expanded];
  CHECK_CUDA_OPTIMIZED(cudaMemcpy(temp_result, ab_mont_gpu,
                                  n_expanded * sizeof(T),
                                  cudaMemcpyDeviceToHost));

  T inv_n_mont = montMod.inv(montMod.from_T(n_expanded));
  for (T i = 0; i < n_expanded; i++) {
    temp_result[i] = montMod.mul(temp_result[i], inv_n_mont);
    temp_result[i] = montMod.to_T(temp_result[i]);
  }

  bit_reverse_permute(temp_result, n_expanded);
  for (T i = 0; i < 2 * n - 1; i++) {
    ab[i] = temp_result[i];
  }

  // 清理
  delete[] a_expanded;
  delete[] b_expanded;
  delete[] a_mont_host;
  delete[] b_mont_host;
  delete[] temp_result;

  cudaStreamDestroy(stream1);
  cudaStreamDestroy(stream2);
}

// =============================================================================
// 优化版本的公开接口
// =============================================================================

template <typename T>
void poly_multiply_ntt_gpu_naive_optimized(T *a, T *b, T *ab, T n, T p,
                                           T omega = 3) {
  poly_multiply_ntt_gpu_optimized<T, MUL_NAIVE>(a, b, ab, n, p, omega);
}

template <typename T>
void poly_multiply_ntt_gpu_mont_optimized(T *a, T *b, T *ab, T n, T p,
                                          T omega = 3) {
  poly_multiply_ntt_gpu_optimized<T, MUL_MONT>(a, b, ab, n, p, omega);
}

template <typename T>
void poly_multiply_ntt_gpu_barrett_optimized(T *a, T *b, T *ab, T n, T p,
                                             T omega = 3) {
  poly_multiply_ntt_gpu_optimized<T, MUL_BARRETT>(a, b, ab, n, p, omega);
}

// 显式实例化
template void poly_multiply_ntt_gpu_naive_optimized<u32>(u32 *a, u32 *b,
                                                         u32 *ab, u32 n, u32 p,
                                                         u32 omega);
template void poly_multiply_ntt_gpu_naive_optimized<u64>(u64 *a, u64 *b,
                                                         u64 *ab, u64 n, u64 p,
                                                         u64 omega);
template void poly_multiply_ntt_gpu_mont_optimized<u32>(u32 *a, u32 *b, u32 *ab,
                                                        u32 n, u32 p,
                                                        u32 omega);
template void poly_multiply_ntt_gpu_mont_optimized<u64>(u64 *a, u64 *b, u64 *ab,
                                                        u64 n, u64 p,
                                                        u64 omega);
template void poly_multiply_ntt_gpu_barrett_optimized<u32>(u32 *a, u32 *b,
                                                           u32 *ab, u32 n,
                                                           u32 p, u32 omega);
template void poly_multiply_ntt_gpu_barrett_optimized<u64>(u64 *a, u64 *b,
                                                           u64 *ab, u64 n,
                                                           u64 p, u64 omega);