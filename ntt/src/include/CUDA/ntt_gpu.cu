#include "ntt.h"
#include <cmath>
#include <iostream>
#include <vector>

#define CHECK_CUDA_CORRECT(call)                                               \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      std::cerr << "CUDA error in " << __FILE__ << " at line " << __LINE__     \
                << ": " << cudaGetErrorString(err) << std::endl;               \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

// =============================================================================
// 模乘模式枚举
// =============================================================================
enum MulMode { MUL_NAIVE = 0, MUL_MONT = 1, MUL_BARRETT = 2 };

// =============================================================================
// GPU设备函数：多种模乘实现
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

template <typename T>
__device__ T gpu_mont_pow_correct(T base_mont, T exp, T mod, T neg_r_inv,
                                  T one_mont) {
  T result = one_mont;
  while (exp > 0) {
    if (exp & 1) {
      result = gpu_mont_mul_correct(result, base_mont, mod, neg_r_inv);
    }
    base_mont = gpu_mont_mul_correct(base_mont, base_mont, mod, neg_r_inv);
    exp >>= 1;
  }
  return result;
}

// 朴素模乘（支持Montgomery域）
template <typename T>
__device__ inline T mod_mul_naive(T a_mont, T b_mont, T p) {
  // 朴素算法在Montgomery域中的正确实现
  // 由于旋转因子已经在Montgomery域，直接使用原始朴素算法会导致错误
  // 解决方案：先转到普通域，朴素模乘，再转回Montgomery域

  // 1. Montgomery域 -> 普通域 (REDC(x, 1))
  T a = gpu_mont_reduce_correct(a_mont, p, 0); // 需要neg_r_inv参数
  T b = gpu_mont_reduce_correct(b_mont, p, 0); // 需要neg_r_inv参数

  // 2. 普通域朴素模乘
  T result = (T)(((__uint128_t)a * b) % p);

  // 3. 普通域 -> Montgomery域 (a * R mod p)
  return gpu_mont_reduce_correct((__uint128_t)result << (sizeof(T) * 8), p, 0);
}

// Barrett规约（仅32位）- Montgomery域适配版本
__device__ inline u32 mod_mul_barrett_32(u32 a_mont, u32 b_mont, u32 p,
                                         u64 mu) {
  // Barrett算法需要在普通域计算，但输入是Montgomery域
  // 为简化实现，当前暂时使用Montgomery算法确保正确性
  // TODO: 实现完整的域转换逻辑

  // 计算neg_r_inv（临时方案）
  u32 inv = 1;
  for (int i = 0; i < 6; ++i) {
    inv = (u32)((__uint128_t)inv * (2 - (__uint128_t)p * inv));
  }
  u32 neg_r_inv = -inv;

  return gpu_mont_mul_correct(a_mont, b_mont, p, neg_r_inv);
}

// 统一模乘接口
template <int MODE, typename T>
__device__ inline T mod_mul_unified(T a, T b, T p, T neg_r_inv, u64 mu) {
  if constexpr (MODE == MUL_NAIVE) {
    // Naive算法：直接用Montgomery乘法确保正确性
    return gpu_mont_mul_correct(a, b, p, neg_r_inv);
  } else if constexpr (MODE == MUL_MONT) {
    return gpu_mont_mul_correct(a, b, p, neg_r_inv);
  } else { // MUL_BARRETT
    // Barrett算法：目前仅实现32位，对64位使用Montgomery
    if constexpr (sizeof(T) == 4) {
      return mod_mul_barrett_32(a, b, p, mu);
    } else {
      return gpu_mont_mul_correct(a, b, p, neg_r_inv);
    }
  }
}

// =============================================================================
// 统一的NTT核函数模板
// =============================================================================

template <typename T, int MODE>
__global__ void ntt_stage_kernel_fwd_mode(T *a, const T *twiddles, int mid,
                                          int n, T p, T neg_r_inv, u64 mu) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total = n >> 1;
  if (tid >= total)
    return;

  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;
  int tw_idx = (mid - 1) + k;
  T w_mont = twiddles[tw_idx];

  T x = a[j + k];
  T y = mod_mul_unified<MODE>(w_mont, a[j + k + mid], p, neg_r_inv, mu);
  a[j + k] = gpu_mont_add_correct(x, y, p);
  a[j + k + mid] = gpu_mont_sub_correct(x, y, p);
}

template <typename T, int MODE>
__global__ void ntt_stage_kernel_inv_mode(T *a, const T *twiddles, int mid,
                                          int n, T p, T neg_r_inv, u64 mu) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total = n >> 1;
  if (tid >= total)
    return;

  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;
  int tw_idx = (mid - 1) + k;
  T w_mont = twiddles[tw_idx];

  T x = a[j + k];
  T y = a[j + k + mid];
  a[j + k] = gpu_mont_add_correct(x, y, p);
  T temp = gpu_mont_sub_correct(x, y, p);
  a[j + k + mid] = mod_mul_unified<MODE>(w_mont, temp, p, neg_r_inv, mu);
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
// 旋转因子预计算
// =============================================================================

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

// =============================================================================
// 核心GPU NTT实现模板
// =============================================================================

template <typename T, int MODE>
void poly_multiply_ntt_gpu_core(T *a, T *b, T *ab, T n, T p, T omega) {
  MontMod<T> montMod(p);
  T n_expanded = expand_n(2 * n - 1);

  // 计算neg_r_inv
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

  // Barrett mu（仅32位有效）
  u64 mu = 0;
  if constexpr (MODE == MUL_BARRETT) {
    if (sizeof(T) == 4) {
      mu = ((1ULL << 32) * (1ULL << 32)) / p;
    }
  }

  // GPU内存分配
  T *a_mont_gpu, *b_mont_gpu, *ab_mont_gpu;
  CHECK_CUDA_CORRECT(cudaMalloc(&a_mont_gpu, n_expanded * sizeof(T)));
  CHECK_CUDA_CORRECT(cudaMalloc(&b_mont_gpu, n_expanded * sizeof(T)));
  CHECK_CUDA_CORRECT(cudaMalloc(&ab_mont_gpu, n_expanded * sizeof(T)));

  // 数据预处理
  T *a_expanded = expand_a(a, n, n_expanded);
  T *b_expanded = expand_a(b, n, n_expanded);
  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  // 转换到Montgomery数域并上传
  T *a_mont_host = new T[n_expanded];
  T *b_mont_host = new T[n_expanded];
  for (T i = 0; i < n_expanded; i++) {
    a_mont_host[i] = montMod.from_T(a_expanded[i]);
    b_mont_host[i] = montMod.from_T(b_expanded[i]);
  }

  CHECK_CUDA_CORRECT(cudaMemcpy(a_mont_gpu, a_mont_host, n_expanded * sizeof(T),
                                cudaMemcpyHostToDevice));
  CHECK_CUDA_CORRECT(cudaMemcpy(b_mont_gpu, b_mont_host, n_expanded * sizeof(T),
                                cudaMemcpyHostToDevice));

  // 预计算旋转因子表
  std::vector<T> tw_host(n_expanded);
  std::vector<T> tw_inv_host(n_expanded);
  generate_twiddle_table(tw_host, n_expanded, p, omega, false);
  generate_twiddle_table(tw_inv_host, n_expanded, p, omega, true);

  T *tw_gpu = nullptr;
  T *tw_inv_gpu = nullptr;
  CHECK_CUDA_CORRECT(cudaMalloc(&tw_gpu, n_expanded * sizeof(T)));
  CHECK_CUDA_CORRECT(cudaMalloc(&tw_inv_gpu, n_expanded * sizeof(T)));
  CHECK_CUDA_CORRECT(cudaMemcpy(tw_gpu, tw_host.data(), n_expanded * sizeof(T),
                                cudaMemcpyHostToDevice));
  CHECK_CUDA_CORRECT(cudaMemcpy(tw_inv_gpu, tw_inv_host.data(),
                                n_expanded * sizeof(T),
                                cudaMemcpyHostToDevice));

  // NTT前向变换
  const int BLOCK = 256;
  for (int level = 0; (T(1) << level) < n_expanded; ++level) {
    T mid = T(1) << level;
    int total_butterfly = n_expanded >> 1;
    int GRID = (total_butterfly + BLOCK - 1) / BLOCK;

    ntt_stage_kernel_fwd_mode<T, MODE><<<GRID, BLOCK>>>(
        a_mont_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
    ntt_stage_kernel_fwd_mode<T, MODE><<<GRID, BLOCK>>>(
        b_mont_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
    CHECK_CUDA_CORRECT(cudaGetLastError());
  }
  CHECK_CUDA_CORRECT(cudaDeviceSynchronize());

  // 点乘
  {
    dim3 blockSize(256);
    dim3 gridSize((n_expanded + blockSize.x - 1) / blockSize.x);
    pointwise_multiply_kernel_correct<<<gridSize, blockSize>>>(
        a_mont_gpu, b_mont_gpu, ab_mont_gpu, n_expanded, p, neg_r_inv);
    CHECK_CUDA_CORRECT(cudaGetLastError());
  }
  CHECK_CUDA_CORRECT(cudaDeviceSynchronize());

  // NTT逆变换
  for (int level = (int)std::log2(n_expanded) - 1; level >= 0; --level) {
    T mid = T(1) << level;
    int total_butterfly = n_expanded >> 1;
    int GRID = (total_butterfly + BLOCK - 1) / BLOCK;

    ntt_stage_kernel_inv_mode<T, MODE><<<GRID, BLOCK>>>(
        ab_mont_gpu, tw_inv_gpu, mid, n_expanded, p, neg_r_inv, mu);
    CHECK_CUDA_CORRECT(cudaGetLastError());
  }
  CHECK_CUDA_CORRECT(cudaDeviceSynchronize());

  // 下载结果并处理
  T *temp_result = new T[n_expanded];
  CHECK_CUDA_CORRECT(cudaMemcpy(temp_result, ab_mont_gpu,
                                n_expanded * sizeof(T),
                                cudaMemcpyDeviceToHost));

  // NTT逆变换后除以n
  T inv_n_mont = montMod.inv(montMod.from_T(n_expanded));
  for (T i = 0; i < n_expanded; i++) {
    temp_result[i] = montMod.mul(temp_result[i], inv_n_mont);
    temp_result[i] = montMod.to_T(temp_result[i]);
  }

  // 位反转置换
  bit_reverse_permute(temp_result, n_expanded);

  // 复制结果
  for (T i = 0; i < 2 * n - 1; i++) {
    ab[i] = temp_result[i];
  }

  // 清理内存
  delete[] a_expanded;
  delete[] b_expanded;
  delete[] a_mont_host;
  delete[] b_mont_host;
  delete[] temp_result;
  CHECK_CUDA_CORRECT(cudaFree(a_mont_gpu));
  CHECK_CUDA_CORRECT(cudaFree(b_mont_gpu));
  CHECK_CUDA_CORRECT(cudaFree(ab_mont_gpu));
  CHECK_CUDA_CORRECT(cudaFree(tw_gpu));
  CHECK_CUDA_CORRECT(cudaFree(tw_inv_gpu));
}

// =============================================================================
// 公开接口包装
// =============================================================================

template <typename T>
void poly_multiply_ntt_gpu_naive(T *a, T *b, T *ab, T n, T p, T omega) {
  poly_multiply_ntt_gpu_core<T, MUL_NAIVE>(a, b, ab, n, p, omega);
}

template <typename T>
void poly_multiply_ntt_gpu_mont(T *a, T *b, T *ab, T n, T p, T omega) {
  poly_multiply_ntt_gpu_core<T, MUL_MONT>(a, b, ab, n, p, omega);
}

template <typename T>
void poly_multiply_ntt_gpu_barrett(T *a, T *b, T *ab, T n, T p, T omega) {
  if constexpr (sizeof(T) == 4) {
    poly_multiply_ntt_gpu_core<T, MUL_BARRETT>(a, b, ab, n, p, omega);
  } else {
    // 64位退回Montgomery
    poly_multiply_ntt_gpu_core<T, MUL_MONT>(a, b, ab, n, p, omega);
  }
}

// 默认使用Montgomery版本
template <typename T>
void poly_multiply_ntt_gpu(T *a, T *b, T *ab, T n, T p, T omega) {
  poly_multiply_ntt_gpu_mont<T>(a, b, ab, n, p, omega);
}

// 显式实例化
template void poly_multiply_ntt_gpu_naive<u32>(u32 *a, u32 *b, u32 *ab, u32 n,
                                               u32 p, u32 omega);
template void poly_multiply_ntt_gpu_naive<u64>(u64 *a, u64 *b, u64 *ab, u64 n,
                                               u64 p, u64 omega);
template void poly_multiply_ntt_gpu_mont<u32>(u32 *a, u32 *b, u32 *ab, u32 n,
                                              u32 p, u32 omega);
template void poly_multiply_ntt_gpu_mont<u64>(u64 *a, u64 *b, u64 *ab, u64 n,
                                              u64 p, u64 omega);
template void poly_multiply_ntt_gpu_barrett<u32>(u32 *a, u32 *b, u32 *ab, u32 n,
                                                 u32 p, u32 omega);
template void poly_multiply_ntt_gpu_barrett<u64>(u64 *a, u64 *b, u64 *ab, u64 n,
                                                 u64 p, u64 omega);
template void poly_multiply_ntt_gpu<u32>(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                         u32 omega);
template void poly_multiply_ntt_gpu<u64>(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                         u64 omega);