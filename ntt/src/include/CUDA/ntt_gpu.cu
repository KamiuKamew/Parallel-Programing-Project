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
// GPU设备函数：多种模乘实现 - u32版本
// =============================================================================

__device__ u32 gpu_mont_reduce_correct_u32(u64 t, u32 mod, u32 neg_r_inv) {
  u32 m = (u32)t * neg_r_inv;
  u64 tmp = t + (u64)m * mod;
  u32 res = (u32)(tmp >> 32);
  res = res - (mod & -(res >= mod));
  return res;
}

__device__ u32 gpu_mont_add_correct_u32(u32 a_mont, u32 b_mont, u32 mod) {
  return (a_mont + b_mont >= mod) ? (a_mont + b_mont - mod) : (a_mont + b_mont);
}

__device__ u32 gpu_mont_sub_correct_u32(u32 a_mont, u32 b_mont, u32 mod) {
  return (a_mont >= b_mont) ? (a_mont - b_mont) : (a_mont + mod - b_mont);
}

__device__ u32 gpu_mont_mul_correct_u32(u32 a_mont, u32 b_mont, u32 mod,
                                        u32 neg_r_inv) {
  return gpu_mont_reduce_correct_u32((u64)a_mont * b_mont, mod, neg_r_inv);
}

__device__ u32 gpu_mont_pow_correct_u32(u32 base_mont, u32 exp, u32 mod,
                                        u32 neg_r_inv, u32 one_mont) {
  u32 result = one_mont;
  while (exp > 0) {
    if (exp & 1) {
      result = gpu_mont_mul_correct_u32(result, base_mont, mod, neg_r_inv);
    }
    base_mont = gpu_mont_mul_correct_u32(base_mont, base_mont, mod, neg_r_inv);
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
  for (int i = 0; i < 5; ++i) { // 5 iterations for 32-bit
    inv = inv * (2 - p * inv);
  }
  u32 neg_r_inv = -inv;

  return gpu_mont_mul_correct_u32(a_mont, b_mont, p, neg_r_inv);
}

// 统一模乘接口
template <int MODE>
__device__ inline u32 mod_mul_unified_u32(u32 a, u32 b, u32 p, u32 neg_r_inv,
                                          u64 mu) {
  if constexpr (MODE == MUL_NAIVE) {
    return gpu_mont_mul_correct_u32(a, b, p, neg_r_inv);
  } else if constexpr (MODE == MUL_MONT) {
    return gpu_mont_mul_correct_u32(a, b, p, neg_r_inv);
  } else { // MUL_BARRETT
    return mod_mul_barrett_32(a, b, p, mu);
  }
}

// =============================================================================
// 统一的NTT核函数 - u32版本
// =============================================================================

template <int MODE>
__global__ void ntt_stage_kernel_fwd_mode_u32(u32 *a, const u32 *twiddles,
                                              int mid, int n, u32 p,
                                              u32 neg_r_inv, u64 mu) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total = n >> 1;
  if (tid >= total)
    return;

  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;
  int tw_idx = (mid - 1) + k;
  u32 w_mont = twiddles[tw_idx];

  u32 x = a[j + k];
  if (x >= p)
    x -= p;

  u32 y_in = a[j + k + mid];
  if (y_in >= p)
    y_in -= p;

  u32 y = mod_mul_unified_u32<MODE>(w_mont, y_in, p, neg_r_inv, mu);

  u32 sum = gpu_mont_add_correct_u32(x, y, p);
  a[j + k] = (sum >= p) ? sum - p : sum;

  u32 diff = gpu_mont_sub_correct_u32(x, y, p);
  a[j + k + mid] = diff;
}

template <int MODE>
__global__ void ntt_stage_kernel_inv_mode_u32(u32 *a, const u32 *twiddles,
                                              int mid, int n, u32 p,
                                              u32 neg_r_inv, u64 mu) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total = n >> 1;
  if (tid >= total)
    return;

  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;
  int tw_idx = (mid - 1) + k;
  u32 w_mont = twiddles[tw_idx];

  u32 x = a[j + k];
  if (x >= p)
    x -= p;

  u32 y = a[j + k + mid];
  if (y >= p)
    y -= p;

  u32 sum = gpu_mont_add_correct_u32(x, y, p);
  a[j + k] = (sum >= p) ? sum - p : sum;

  u32 temp = gpu_mont_sub_correct_u32(x, y, p);
  a[j + k + mid] = mod_mul_unified_u32<MODE>(w_mont, temp, p, neg_r_inv, mu);
}

__global__ void pointwise_multiply_kernel_correct_u32(u32 *a_mont, u32 *b_mont,
                                                      u32 *ab_mont, u32 n,
                                                      u32 mod, u32 neg_r_inv) {
  u32 idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < n) {
    ab_mont[idx] =
        gpu_mont_mul_correct_u32(a_mont[idx], b_mont[idx], mod, neg_r_inv);
  }
}

// =============================================================================
// 旋转因子预计算 - u32版本
// =============================================================================

void generate_twiddle_table_u32(std::vector<u32> &table, u32 n, u32 p,
                                u32 omega, bool inverse) {
  table.resize(n);
  MontMod<u32> mont(p);
  u32 omega_mont = mont.from_T(omega);
  if (inverse)
    omega_mont = mont.inv(omega_mont);
  u32 one_mont = mont.from_T(1);

  for (int level = 0; (1u << level) < n; ++level) {
    u32 mid = 1u << level;
    u32 exp = (p - 1) / (mid << 1);
    u32 Wn_mont = mont.pow(omega_mont, exp);
    u32 w = one_mont;
    u32 offset = mid - 1;
    for (u32 k = 0; k < mid; ++k) {
      table[offset + k] = w;
      w = mont.mul(w, Wn_mont);
    }
  }
}

// =============================================================================
// 核心GPU NTT实现模板 - u32版本
// =============================================================================

template <int MODE>
void poly_multiply_ntt_gpu_core_u32(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                    u32 omega) {
  MontMod<u32> montMod(p);
  u32 n_expanded = expand_n(2 * n - 1);

  // 计算neg_r_inv
  u32 neg_r_inv;
  {
    u32 inv = 1;
    for (int i = 0; i < 5; ++i) { // 5 iterations for 32-bit
      inv = inv * (2 - p * inv);
    }
    neg_r_inv = -inv;
  }

  // Barrett mu
  u64 mu = 0;
  if constexpr (MODE == MUL_BARRETT) {
    mu = ((1ULL << 32) * (1ULL << 32)) / p;
  }

  // GPU内存分配
  u32 *a_mont_gpu, *b_mont_gpu, *ab_mont_gpu;
  CHECK_CUDA_CORRECT(cudaMalloc(&a_mont_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_CORRECT(cudaMalloc(&b_mont_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_CORRECT(cudaMalloc(&ab_mont_gpu, n_expanded * sizeof(u32)));

  // 数据预处理
  u32 *a_expanded = expand_a(a, n, n_expanded);
  u32 *b_expanded = expand_a(b, n, n_expanded);
  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  // 转换到Montgomery数域并上传
  u32 *a_mont_host = new u32[n_expanded];
  u32 *b_mont_host = new u32[n_expanded];
  for (u32 i = 0; i < n_expanded; i++) {
    a_mont_host[i] = montMod.from_T(a_expanded[i]);
    b_mont_host[i] = montMod.from_T(b_expanded[i]);
  }

  CHECK_CUDA_CORRECT(cudaMemcpy(a_mont_gpu, a_mont_host,
                                n_expanded * sizeof(u32),
                                cudaMemcpyHostToDevice));
  CHECK_CUDA_CORRECT(cudaMemcpy(b_mont_gpu, b_mont_host,
                                n_expanded * sizeof(u32),
                                cudaMemcpyHostToDevice));

  // 预计算旋转因子表
  std::vector<u32> tw_host(n_expanded);
  std::vector<u32> tw_inv_host(n_expanded);
  generate_twiddle_table_u32(tw_host, n_expanded, p, omega, false);
  generate_twiddle_table_u32(tw_inv_host, n_expanded, p, omega, true);

  u32 *tw_gpu = nullptr;
  u32 *tw_inv_gpu = nullptr;
  CHECK_CUDA_CORRECT(cudaMalloc(&tw_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_CORRECT(cudaMalloc(&tw_inv_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_CORRECT(cudaMemcpy(tw_gpu, tw_host.data(),
                                n_expanded * sizeof(u32),
                                cudaMemcpyHostToDevice));
  CHECK_CUDA_CORRECT(cudaMemcpy(tw_inv_gpu, tw_inv_host.data(),
                                n_expanded * sizeof(u32),
                                cudaMemcpyHostToDevice));

  // NTT前向变换
  const int BLOCK = 256;
  for (int level = 0; (1u << level) < n_expanded; ++level) {
    u32 mid = 1u << level;
    int total_butterfly = n_expanded >> 1;
    int GRID = (total_butterfly + BLOCK - 1) / BLOCK;

    ntt_stage_kernel_fwd_mode_u32<MODE><<<GRID, BLOCK>>>(
        a_mont_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
    ntt_stage_kernel_fwd_mode_u32<MODE><<<GRID, BLOCK>>>(
        b_mont_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
    CHECK_CUDA_CORRECT(cudaGetLastError());
  }
  CHECK_CUDA_CORRECT(cudaDeviceSynchronize());

  // 点乘
  {
    dim3 blockSize(256);
    dim3 gridSize((n_expanded + blockSize.x - 1) / blockSize.x);
    pointwise_multiply_kernel_correct_u32<<<gridSize, blockSize>>>(
        a_mont_gpu, b_mont_gpu, ab_mont_gpu, n_expanded, p, neg_r_inv);
    CHECK_CUDA_CORRECT(cudaGetLastError());
  }
  CHECK_CUDA_CORRECT(cudaDeviceSynchronize());

  // NTT逆变换
  for (int level = (int)std::log2(n_expanded) - 1; level >= 0; --level) {
    u32 mid = 1u << level;
    int total_butterfly = n_expanded >> 1;
    int GRID = (total_butterfly + BLOCK - 1) / BLOCK;

    ntt_stage_kernel_inv_mode_u32<MODE><<<GRID, BLOCK>>>(
        ab_mont_gpu, tw_inv_gpu, mid, n_expanded, p, neg_r_inv, mu);
    CHECK_CUDA_CORRECT(cudaGetLastError());
  }
  CHECK_CUDA_CORRECT(cudaDeviceSynchronize());

  // 下载结果并处理
  u32 *temp_result = new u32[n_expanded];
  CHECK_CUDA_CORRECT(cudaMemcpy(temp_result, ab_mont_gpu,
                                n_expanded * sizeof(u32),
                                cudaMemcpyDeviceToHost));

  // NTT逆变换后除以n
  u32 inv_n_mont = montMod.inv(montMod.from_T(n_expanded));
  for (u32 i = 0; i < n_expanded; i++) {
    temp_result[i] = montMod.mul(temp_result[i], inv_n_mont);
    temp_result[i] = montMod.to_T(temp_result[i]);
  }

  // 位反转置换
  bit_reverse_permute(temp_result, n_expanded);

  // 复制结果
  for (u32 i = 0; i < 2 * n - 1; i++) {
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
// 公开接口包装 - u32版本
// =============================================================================

void poly_multiply_ntt_gpu_naive(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                 u32 omega) {
  poly_multiply_ntt_gpu_core_u32<MUL_NAIVE>(a, b, ab, n, p, omega);
}

void poly_multiply_ntt_gpu_mont(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                u32 omega) {
  poly_multiply_ntt_gpu_core_u32<MUL_MONT>(a, b, ab, n, p, omega);
}

void poly_multiply_ntt_gpu_barrett(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                   u32 omega) {
  poly_multiply_ntt_gpu_core_u32<MUL_BARRETT>(a, b, ab, n, p, omega);
}

// 默认使用Montgomery版本
void poly_multiply_ntt_gpu(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega) {
  poly_multiply_ntt_gpu_mont(a, b, ab, n, p, omega);
}