#include "ntt.h"
#include <cmath>
#include <iostream>
#include <vector>

// CUDA错误检查宏
#define CHECK_CUDA_STAGE2(call)                                                \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,         \
              cudaGetErrorString(err));                                        \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)

// 模乘模式枚举
enum Stage2MulMode {
  STAGE2_MUL_NAIVE = 0,
  STAGE2_MUL_MONT = 1,
  STAGE2_MUL_BARRETT = 2
};

// =============================================================================
// 旋转因子预计算函数
// =============================================================================
void generate_twiddle_table_stage2(std::vector<u32> &table, u32 n, u32 p,
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
// 优化的GPU设备函数 - Stage 2
// =============================================================================
__device__ u32 gpu_stage2_mont_reduce(u64 t, u32 mod, u32 neg_r_inv) {
  u32 m = (u32)t * neg_r_inv;
  u64 tmp = t + (u64)m * mod;
  u32 res = (u32)(tmp >> 32);
  return res - (mod & -(res >= mod));
}

__device__ u32 gpu_stage2_mont_mul(u32 a, u32 b, u32 mod, u32 neg_r_inv) {
  return gpu_stage2_mont_reduce((u64)a * b, mod, neg_r_inv);
}

__device__ u32 gpu_stage2_mont_add(u32 a, u32 b, u32 mod) {
  return (a + b >= mod) ? (a + b - mod) : (a + b);
}

__device__ u32 gpu_stage2_mont_sub(u32 a, u32 b, u32 mod) {
  return (a >= b) ? (a - b) : (a + mod - b);
}

// Barrett规约优化版本
__device__ inline u32 mod_mul_barrett_stage2(u32 a, u32 b, u32 p, u64 mu) {
  u64 product = (u64)a * b;
  u64 q = (product * mu) >> 32;
  u32 r = (u32)(product - q * p);
  return r >= p ? r - p : r;
}

// 统一模乘接口 - Stage 2优化
template <int MODE>
__device__ inline u32 mod_mul_unified_stage2(u32 a, u32 b, u32 p, u32 neg_r_inv,
                                             u64 mu) {
  if constexpr (MODE == STAGE2_MUL_NAIVE) {
    return ((u64)a * b) % p;
  } else if constexpr (MODE == STAGE2_MUL_MONT) {
    return gpu_stage2_mont_mul(a, b, p, neg_r_inv);
  } else { // STAGE2_MUL_BARRETT
    return mod_mul_barrett_stage2(a, b, p, mu);
  }
}

// =============================================================================
// Stage 2优化核函数 - 标准NTT阶段（未融合）
// =============================================================================
template <int MODE>
__global__ void ntt_stage2_kernel_fwd(u32 *a, const u32 *twiddles, int mid,
                                      int n, u32 p, u32 neg_r_inv, u64 mu) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total_butterflies = n >> 1;

  if (tid >= total_butterflies)
    return;

  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;

  int tw_idx = (mid - 1) + k;
  int idx1 = j + k;
  int idx2 = idx1 + mid;

  u32 w_mont = twiddles[tw_idx];
  u32 x = a[idx1];
  u32 y_data = a[idx2];

  u32 y = mod_mul_unified_stage2<MODE>(w_mont, y_data, p, neg_r_inv, mu);

  u32 sum = gpu_stage2_mont_add(x, y, p);
  u32 diff = gpu_stage2_mont_sub(x, y, p);

  a[idx1] = sum;
  a[idx2] = diff;
}

// =============================================================================
// Stage 2核心优化：融合的最后NTT阶段 + 点乘核函数
// =============================================================================
template <int MODE>
__global__ void ntt_fused_final_stage_and_pointwise(u32 *a, u32 *b, u32 *ab,
                                                    const u32 *twiddles,
                                                    int final_mid, int n, u32 p,
                                                    u32 neg_r_inv, u64 mu) {

  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total_butterflies = n >> 1;

  if (tid >= total_butterflies)
    return;

  // 执行最后一个NTT阶段的蝶形操作
  int step = final_mid << 1;
  int j = (tid / final_mid) * step;
  int k = tid % final_mid;

  int tw_idx = (final_mid - 1) + k;
  int idx1 = j + k;
  int idx2 = idx1 + final_mid;

  u32 w_mont = twiddles[tw_idx];

  // 对a数组执行最后的NTT阶段
  u32 a_x = a[idx1];
  u32 a_y_data = a[idx2];
  u32 a_y = mod_mul_unified_stage2<MODE>(w_mont, a_y_data, p, neg_r_inv, mu);
  u32 a_sum = gpu_stage2_mont_add(a_x, a_y, p);
  u32 a_diff = gpu_stage2_mont_sub(a_x, a_y, p);

  // 对b数组执行最后的NTT阶段
  u32 b_x = b[idx1];
  u32 b_y_data = b[idx2];
  u32 b_y = mod_mul_unified_stage2<MODE>(w_mont, b_y_data, p, neg_r_inv, mu);
  u32 b_sum = gpu_stage2_mont_add(b_x, b_y, p);
  u32 b_diff = gpu_stage2_mont_sub(b_x, b_y, p);

  // Stage 2核心优化：立即执行点乘，避免写回和重新读取
  ab[idx1] = gpu_stage2_mont_mul(a_sum, b_sum, p, neg_r_inv);
  ab[idx2] = gpu_stage2_mont_mul(a_diff, b_diff, p, neg_r_inv);
}

// 标准逆NTT核函数
template <int MODE>
__global__ void ntt_stage2_kernel_inv(u32 *a, const u32 *twiddles_inv, int mid,
                                      int n, u32 p, u32 neg_r_inv, u64 mu) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total_butterflies = n >> 1;

  if (tid >= total_butterflies)
    return;

  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;

  int tw_idx = (mid - 1) + k;
  int idx1 = j + k;
  int idx2 = idx1 + mid;

  u32 w_mont = twiddles_inv[tw_idx];
  u32 x = a[idx1];
  u32 y = a[idx2];

  u32 sum = gpu_stage2_mont_add(x, y, p);
  u32 temp = gpu_stage2_mont_sub(x, y, p);
  u32 diff = mod_mul_unified_stage2<MODE>(w_mont, temp, p, neg_r_inv, mu);

  a[idx1] = sum;
  a[idx2] = diff;
}

// =============================================================================
// Stage 2优化主函数
// =============================================================================
template <int MODE>
void poly_multiply_ntt_gpu_stage2_core(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                       u32 omega) {
  MontMod<u32> montMod(p);
  u32 n_expanded = expand_n(2 * n - 1);

  // 计算Montgomery参数
  u32 neg_r_inv;
  {
    u32 inv = 1;
    for (int i = 0; i < 5; ++i) {
      inv = inv * (2 - p * inv);
    }
    neg_r_inv = -inv;
  }

  // Barrett参数
  u64 mu = 0;
  if constexpr (MODE == STAGE2_MUL_BARRETT) {
    mu = ((unsigned __int128)1 << 64) / p;
  }

  // 预计算旋转因子表
  std::vector<u32> tw_host(n_expanded);
  std::vector<u32> tw_inv_host(n_expanded);
  generate_twiddle_table_stage2(tw_host, n_expanded, p, omega, false);
  generate_twiddle_table_stage2(tw_inv_host, n_expanded, p, omega, true);

  // 上传旋转因子到GPU
  u32 *tw_gpu, *tw_inv_gpu;
  CHECK_CUDA_STAGE2(cudaMalloc(&tw_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_STAGE2(cudaMalloc(&tw_inv_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_STAGE2(cudaMemcpy(tw_gpu, tw_host.data(), n_expanded * sizeof(u32),
                               cudaMemcpyHostToDevice));
  CHECK_CUDA_STAGE2(cudaMemcpy(tw_inv_gpu, tw_inv_host.data(),
                               n_expanded * sizeof(u32),
                               cudaMemcpyHostToDevice));

  // GPU内存分配
  u32 *a_gpu, *b_gpu, *ab_gpu;
  CHECK_CUDA_STAGE2(cudaMalloc(&a_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_STAGE2(cudaMalloc(&b_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_STAGE2(cudaMalloc(&ab_gpu, n_expanded * sizeof(u32)));

  // 数据预处理
  u32 *a_expanded = expand_a(a, n, n_expanded);
  u32 *b_expanded = expand_a(b, n, n_expanded);
  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  // 转换到Montgomery域
  u32 *a_mont_host = new u32[n_expanded];
  u32 *b_mont_host = new u32[n_expanded];
  for (u32 i = 0; i < n_expanded; i++) {
    a_mont_host[i] = montMod.from_T(a_expanded[i]);
    b_mont_host[i] = montMod.from_T(b_expanded[i]);
  }

  // 上传数据到GPU
  CHECK_CUDA_STAGE2(cudaMemcpy(a_gpu, a_mont_host, n_expanded * sizeof(u32),
                               cudaMemcpyHostToDevice));
  CHECK_CUDA_STAGE2(cudaMemcpy(b_gpu, b_mont_host, n_expanded * sizeof(u32),
                               cudaMemcpyHostToDevice));

  // Stage 2优化：使用256线程块
  const int BLOCK_SIZE = 256;

  // 前向NTT变换 - 除了最后一个阶段
  u32 final_mid = n_expanded >> 1; // 最后一个mid值
  for (u32 mid = 1; mid < final_mid; mid <<= 1) {
    int total_butterflies = n_expanded >> 1;
    int grid_size = (total_butterflies + BLOCK_SIZE - 1) / BLOCK_SIZE;

    ntt_stage2_kernel_fwd<MODE><<<grid_size, BLOCK_SIZE>>>(
        a_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
    ntt_stage2_kernel_fwd<MODE><<<grid_size, BLOCK_SIZE>>>(
        b_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);

    CHECK_CUDA_STAGE2(cudaGetLastError());
  }
  CHECK_CUDA_STAGE2(cudaDeviceSynchronize());

  // Stage 2核心优化：融合最后的NTT阶段和点乘操作
  {
    int total_butterflies = n_expanded >> 1;
    int grid_size = (total_butterflies + BLOCK_SIZE - 1) / BLOCK_SIZE;

    ntt_fused_final_stage_and_pointwise<MODE><<<grid_size, BLOCK_SIZE>>>(
        a_gpu, b_gpu, ab_gpu, tw_gpu, final_mid, n_expanded, p, neg_r_inv, mu);

    CHECK_CUDA_STAGE2(cudaGetLastError());
  }
  CHECK_CUDA_STAGE2(cudaDeviceSynchronize());

  // 逆NTT变换
  for (u32 mid = n_expanded >> 1; mid >= 1; mid >>= 1) {
    int total_butterflies = n_expanded >> 1;
    int grid_size = (total_butterflies + BLOCK_SIZE - 1) / BLOCK_SIZE;

    ntt_stage2_kernel_inv<MODE><<<grid_size, BLOCK_SIZE>>>(
        ab_gpu, tw_inv_gpu, mid, n_expanded, p, neg_r_inv, mu);

    CHECK_CUDA_STAGE2(cudaGetLastError());
  }
  CHECK_CUDA_STAGE2(cudaDeviceSynchronize());

  // 下载结果
  u32 *ab_mont_host = new u32[n_expanded];
  CHECK_CUDA_STAGE2(cudaMemcpy(ab_mont_host, ab_gpu, n_expanded * sizeof(u32),
                               cudaMemcpyDeviceToHost));

  // 位反转置换
  bit_reverse_permute(ab_mont_host, n_expanded);

  // 逆变换后处理
  u32 inv_n_mont = montMod.inv(montMod.from_T(n_expanded));
  for (u32 i = 0; i < 2 * n - 1; i++) {
    ab[i] = montMod.to_T(montMod.mul(ab_mont_host[i], inv_n_mont));
  }

  // 清理内存
  delete[] a_expanded;
  delete[] b_expanded;
  delete[] a_mont_host;
  delete[] b_mont_host;
  delete[] ab_mont_host;
  CHECK_CUDA_STAGE2(cudaFree(a_gpu));
  CHECK_CUDA_STAGE2(cudaFree(b_gpu));
  CHECK_CUDA_STAGE2(cudaFree(ab_gpu));
  CHECK_CUDA_STAGE2(cudaFree(tw_gpu));
  CHECK_CUDA_STAGE2(cudaFree(tw_inv_gpu));
}

// =============================================================================
// 公开接口 - Stage 2优化版本
// =============================================================================
void poly_multiply_ntt_gpu_stage2_mont(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                       u32 omega) {
  poly_multiply_ntt_gpu_stage2_core<STAGE2_MUL_MONT>(a, b, ab, n, p, omega);
}

void poly_multiply_ntt_gpu_stage2_naive(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                        u32 omega) {
  poly_multiply_ntt_gpu_stage2_core<STAGE2_MUL_NAIVE>(a, b, ab, n, p, omega);
}

void poly_multiply_ntt_gpu_stage2_barrett(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                          u32 omega) {
  poly_multiply_ntt_gpu_stage2_core<STAGE2_MUL_BARRETT>(a, b, ab, n, p, omega);
}