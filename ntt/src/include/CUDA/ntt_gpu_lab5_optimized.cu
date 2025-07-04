#include "ntt.h"
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// CUDA错误检查宏
#define CHECK_CUDA_OPTIMIZED(call)                                             \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,         \
              cudaGetErrorString(err));                                        \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)

// 模乘模式枚举 - 保持与ntt_gpu.cu一致
enum MulMode { OPT_MUL_NAIVE = 0, OPT_MUL_MONT = 1, OPT_MUL_BARRETT = 2 };

// =============================================================================
// 统一的GPU设备函数 - u32版本
// =============================================================================
__device__ u32 gpu_opt_mont_reduce(u64 t, u32 mod, u32 neg_r_inv) {
  u32 m = (u32)t * neg_r_inv;
  u64 tmp = t + (u64)m * mod;
  u32 res = (u32)(tmp >> 32);
  res = res - (mod & -(res >= mod));
  return res;
}

__device__ u32 gpu_opt_mont_mul(u32 a, u32 b, u32 mod, u32 neg_r_inv) {
  return gpu_opt_mont_reduce((u64)a * b, mod, neg_r_inv);
}

__device__ u32 gpu_opt_mont_add(u32 a, u32 b, u32 mod) {
  return (a + b >= mod) ? (a + b - mod) : (a + b);
}

__device__ u32 gpu_opt_mont_sub(u32 a, u32 b, u32 mod) {
  return (a >= b) ? (a - b) : (a + mod - b);
}

__device__ inline u32 mod_mul_barrett_opt_32(u32 a, u32 b, u32 p, u64 mu) {
  u32 t = (u32)(((unsigned __int128)a * b) >> 32);
  u32 r = a * b - t * p;
  return r >= p ? r - p : r;
}

template <int MODE>
__device__ inline u32 mod_mul_unified_opt(u32 a, u32 b, u32 p, u32 neg_r_inv,
                                          u64 mu) {
  if constexpr (MODE == OPT_MUL_NAIVE) {
    return ((u64)a * b) % p;
  } else if constexpr (MODE == OPT_MUL_MONT) {
    return gpu_opt_mont_mul(a, b, p, neg_r_inv);
  } else { // OPT_MUL_BARRETT
    return mod_mul_barrett_opt_32(a, b, p, mu);
  }
}

// =============================================================================
// 优化的NTT核函数 - u32版本
// =============================================================================

// 共享内存符号
extern __shared__ u32 shared_twiddles[];

template <int MODE>
__global__ void ntt_stage_kernel_fwd_optimized_u32(u32 *a, const u32 *twiddles,
                                                   int mid, int n, u32 p,
                                                   u32 neg_r_inv, u64 mu) {
  int tid_global = blockIdx.x * blockDim.x + threadIdx.x;
  int total_butterflies = n >> 1;
  if (tid_global >= total_butterflies)
    return;

  // Load twiddles to shared memory
  if (mid <= blockDim.x) {
    if (threadIdx.x < mid) {
      shared_twiddles[threadIdx.x] = twiddles[(mid - 1) + threadIdx.x];
    }
    __syncthreads();
  }

  int step = mid << 1;
  int j = (tid_global / mid) * step;
  int k = tid_global % mid;

  u32 w_mont;
  if (mid <= blockDim.x) {
    w_mont = shared_twiddles[k];
  } else {
    int tw_idx = (mid - 1) + k;
    w_mont = twiddles[tw_idx];
  }

  int index1 = j + k;
  int index2 = index1 + mid;

  u32 x = a[index1];
  if (x >= p)
    x -= p;

  u32 y_in = a[index2];
  if (y_in >= p)
    y_in -= p;

  u32 y = mod_mul_unified_opt<MODE>(w_mont, y_in, p, neg_r_inv, mu);

  u32 sum = gpu_opt_mont_add(x, y, p);
  a[index1] = (sum >= p) ? sum - p : sum;

  u32 diff = gpu_opt_mont_sub(x, y, p);
  a[index2] = diff;
}

template <int MODE>
__global__ void ntt_stage_kernel_inv_optimized_u32(u32 *a, const u32 *twiddles,
                                                   int mid, int n, u32 p,
                                                   u32 neg_r_inv, u64 mu) {
  int tid_global = blockIdx.x * blockDim.x + threadIdx.x;
  int total_butterflies = n >> 1;
  if (tid_global >= total_butterflies)
    return;

  // Load twiddles to shared memory
  if (mid <= blockDim.x) {
    if (threadIdx.x < mid) {
      shared_twiddles[threadIdx.x] = twiddles[(mid - 1) + threadIdx.x];
    }
    __syncthreads();
  }

  int step = mid << 1;
  int j = (tid_global / mid) * step;
  int k = tid_global % mid;

  u32 w_mont;
  if (mid <= blockDim.x) {
    w_mont = shared_twiddles[k];
  } else {
    int tw_idx = (mid - 1) + k;
    w_mont = twiddles[tw_idx];
  }

  int index1 = j + k;
  int index2 = index1 + mid;

  u32 x = a[index1];
  if (x >= p)
    x -= p;

  u32 y = a[index2];
  if (y >= p)
    y -= p;

  u32 sum = gpu_opt_mont_add(x, y, p);
  a[index1] = (sum >= p) ? sum - p : sum;

  u32 temp = gpu_opt_mont_sub(x, y, p);
  a[index2] = mod_mul_unified_opt<MODE>(w_mont, temp, p, neg_r_inv, mu);
}

__global__ void pointwise_multiply_kernel_optimized_u32(u32 *a, u32 *b, u32 *ab,
                                                        u32 n, u32 p,
                                                        u32 neg_r_inv) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    ab[idx] = gpu_opt_mont_mul(a[idx], b[idx], p, neg_r_inv);
  }
}

// =============================================================================
// 主机端核心实现 - u32版本
// =============================================================================

template <int MODE>
void poly_multiply_ntt_gpu_optimized_core(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                          u32 omega) {
  MontMod<u32> montMod(p);
  u32 n_expanded = expand_n(2 * n - 1);

  u32 neg_r_inv;
  {
    u32 inv = 1;
    for (int i = 0; i < 5; ++i) { // 5 iterations for 32-bit
      inv = inv * (2 - p * inv);
    }
    neg_r_inv = -inv;
  }

  u64 mu = 0;
  if constexpr (MODE == OPT_MUL_BARRETT) {
    mu = ((unsigned __int128)1 << 64) / p;
  }

  u32 *a_gpu, *b_gpu, *ab_gpu;
  CHECK_CUDA_OPTIMIZED(cudaMalloc(&a_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_OPTIMIZED(cudaMalloc(&b_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_OPTIMIZED(cudaMalloc(&ab_gpu, n_expanded * sizeof(u32)));

  u32 *a_expanded = expand_a(a, n, n_expanded);
  u32 *b_expanded = expand_a(b, n, n_expanded);

  u32 *a_mont_host = new u32[n_expanded];
  u32 *b_mont_host = new u32[n_expanded];
  for (u32 i = 0; i < n_expanded; i++) {
    a_mont_host[i] = montMod.from_T(a_expanded[i]);
    b_mont_host[i] = montMod.from_T(b_expanded[i]);
  }

  CHECK_CUDA_OPTIMIZED(cudaMemcpy(a_gpu, a_mont_host, n_expanded * sizeof(u32),
                                  cudaMemcpyHostToDevice));
  CHECK_CUDA_OPTIMIZED(cudaMemcpy(b_gpu, b_mont_host, n_expanded * sizeof(u32),
                                  cudaMemcpyHostToDevice));

  std::vector<u32> tw_host(n_expanded), tw_inv_host(n_expanded);
  generate_twiddle_table_u32(tw_host, n_expanded, p, omega, false);
  generate_twiddle_table_u32(tw_inv_host, n_expanded, p, omega, true);

  u32 *tw_gpu, *tw_inv_gpu;
  CHECK_CUDA_OPTIMIZED(cudaMalloc(&tw_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_OPTIMIZED(cudaMalloc(&tw_inv_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_OPTIMIZED(cudaMemcpy(tw_gpu, tw_host.data(),
                                  n_expanded * sizeof(u32),
                                  cudaMemcpyHostToDevice));
  CHECK_CUDA_OPTIMIZED(cudaMemcpy(tw_inv_gpu, tw_inv_host.data(),
                                  n_expanded * sizeof(u32),
                                  cudaMemcpyHostToDevice));

  // 自适应线程配置
  int block_size = (n_expanded <= 1024) ? 128 : 256;
  if (n_expanded > 65536)
    block_size = 512;

  bit_reverse_permute(a_mont_host, n_expanded);
  bit_reverse_permute(b_mont_host, n_expanded);
  CHECK_CUDA_OPTIMIZED(cudaMemcpy(a_gpu, a_mont_host, n_expanded * sizeof(u32),
                                  cudaMemcpyHostToDevice));
  CHECK_CUDA_OPTIMIZED(cudaMemcpy(b_gpu, b_mont_host, n_expanded * sizeof(u32),
                                  cudaMemcpyHostToDevice));

  for (u32 mid = 1; mid < n_expanded; mid <<= 1) {
    int total_butterflies = n_expanded >> 1;
    int grid_size = (total_butterflies + block_size - 1) / block_size;
    size_t shared_mem_size = (mid <= block_size) ? mid * sizeof(u32) : 0;
    ntt_stage_kernel_fwd_optimized_u32<MODE>
        <<<grid_size, block_size, shared_mem_size>>>(
            a_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
    ntt_stage_kernel_fwd_optimized_u32<MODE>
        <<<grid_size, block_size, shared_mem_size>>>(
            b_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
  }
  CHECK_CUDA_OPTIMIZED(cudaDeviceSynchronize());

  int grid_size_pm = (n_expanded + block_size - 1) / block_size;
  pointwise_multiply_kernel_optimized_u32<<<grid_size_pm, block_size>>>(
      a_gpu, b_gpu, ab_gpu, n_expanded, p, neg_r_inv);
  CHECK_CUDA_OPTIMIZED(cudaDeviceSynchronize());

  for (u32 mid = n_expanded >> 1; mid >= 1; mid >>= 1) {
    int total_butterflies = n_expanded >> 1;
    int grid_size = (total_butterflies + block_size - 1) / block_size;
    size_t shared_mem_size = (mid <= block_size) ? mid * sizeof(u32) : 0;
    ntt_stage_kernel_inv_optimized_u32<MODE>
        <<<grid_size, block_size, shared_mem_size>>>(
            ab_gpu, tw_inv_gpu, mid, n_expanded, p, neg_r_inv, mu);
  }
  CHECK_CUDA_OPTIMIZED(cudaDeviceSynchronize());

  u32 *ab_mont_host = new u32[n_expanded];
  CHECK_CUDA_OPTIMIZED(cudaMemcpy(
      ab_mont_host, ab_gpu, n_expanded * sizeof(u32), cudaMemcpyDeviceToHost));

  bit_reverse_permute(ab_mont_host, n_expanded);

  u32 inv_n_mont = montMod.inv(montMod.from_T(n_expanded));
  for (u32 i = 0; i < 2 * n - 1; i++) {
    ab[i] = montMod.to_T(montMod.mul(ab_mont_host[i], inv_n_mont));
  }

  delete[] a_expanded;
  delete[] b_expanded;
  delete[] a_mont_host;
  delete[] b_mont_host;
  delete[] ab_mont_host;
  CHECK_CUDA_OPTIMIZED(cudaFree(a_gpu));
  CHECK_CUDA_OPTIMIZED(cudaFree(b_gpu));
  CHECK_CUDA_OPTIMIZED(cudaFree(ab_gpu));
  CHECK_CUDA_OPTIMIZED(cudaFree(tw_gpu));
  CHECK_CUDA_OPTIMIZED(cudaFree(tw_inv_gpu));
}

// =============================================================================
// 公开接口 - u32版本
// =============================================================================
void poly_multiply_ntt_gpu_mont_optimized(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                          u32 omega) {
  poly_multiply_ntt_gpu_optimized_core<OPT_MUL_MONT>(a, b, ab, n, p, omega);
}

void poly_multiply_ntt_gpu_naive_optimized(u32 *a, u32 *b, u32 *ab, u32 n,
                                           u32 p, u32 omega) {
  poly_multiply_ntt_gpu_optimized_core<OPT_MUL_NAIVE>(a, b, ab, n, p, omega);
}

void poly_multiply_ntt_gpu_barrett_optimized(u32 *a, u32 *b, u32 *ab, u32 n,
                                             u32 p, u32 omega) {
  poly_multiply_ntt_gpu_optimized_core<OPT_MUL_BARRETT>(a, b, ab, n, p, omega);
}