#include "ntt.h"
#include <cmath>
#include <iostream>
#include <vector>

// CUDA错误检查宏
#define CHECK_CUDA_STAGE3(call)                                                \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,         \
              cudaGetErrorString(err));                                        \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)

// 模乘模式枚举
enum Stage3MulMode {
  STAGE1_MUL_NAIVE = 0,
  STAGE1_MUL_MONT = 1,
  STAGE1_MUL_BARRETT = 2
};

// =============================================================================
// 旋转因子预计算函数
// =============================================================================
void generate_twiddle_table_stage3(std::vector<u32> &table, u32 n, u32 p,
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
// 优化的GPU设备函数 - Stage 1
// =============================================================================
__device__ __forceinline__ u32 gpu_stage3_mont_reduce(u64 t, u32 mod,
                                                      u32 neg_r_inv) {
  u32 m = (u32)t * neg_r_inv;
  u64 tmp = t + (u64)m * mod;
  u32 res = (u32)(tmp >> 32);
  return res - (mod & -(res >= mod));
}

__device__ __forceinline__ u32 gpu_stage3_mont_mul(u32 a, u32 b, u32 mod,
                                                   u32 neg_r_inv) {
  return gpu_stage3_mont_reduce((u64)a * b, mod, neg_r_inv);
}

__device__ __forceinline__ u32 gpu_stage3_mont_add(u32 a, u32 b, u32 mod) {
  return (a + b >= mod) ? (a + b - mod) : (a + b);
}

__device__ __forceinline__ u32 gpu_stage3_mont_sub(u32 a, u32 b, u32 mod) {
  return (a >= b) ? (a - b) : (a + mod - b);
}

// Barrett规约优化版本
__device__ inline u32 mod_mul_barrett_stage3(u32 a, u32 b, u32 p, u64 mu) {
  u64 product = (u64)a * b;
  u64 q = (product * mu) >> 32;
  u32 r = (u32)(product - q * p);
  return r >= p ? r - p : r;
}

// 统一模乘接口 - Stage 1优化
template <int MODE>
__device__ inline u32 mod_mul_unified_stage3(u32 a, u32 b, u32 p, u32 neg_r_inv,
                                             u64 mu) {
  if constexpr (MODE == STAGE1_MUL_NAIVE) {
    return ((u64)a * b) % p;
  } else if constexpr (MODE == STAGE1_MUL_MONT) {
    return gpu_stage3_mont_mul(a, b, p, neg_r_inv);
  } else { // STAGE1_MUL_BARRETT
    return mod_mul_barrett_stage3(a, b, p, mu);
  }
}

// =============================================================================
// Stage 1优化核函数 - 优化内存访问模式 + 寄存器重用
// =============================================================================
template <int MODE>
__global__ void ntt_stage3_kernel_fwd(u32 *a, const u32 *twiddles, int mid,
                                      int n, u32 p, u32 neg_r_inv, u64 mu) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total_butterflies = n >> 1;

  if (tid >= total_butterflies)
    return;

  // 计算蝶形操作的位置
  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;

  // Stage 1优化：预计算索引，减少重复计算
  int tw_idx = (mid - 1) + k;
  int idx1 = j + k;
  int idx2 = idx1 + mid;

  // 寄存器重用优化：一次性加载所有需要的数据
  u32 w_mont = twiddles[tw_idx];
  u32 x = a[idx1];
  u32 y_data = a[idx2];

  // 执行蝶形操作 - 使用寄存器变量
  u32 y = mod_mul_unified_stage3<MODE>(w_mont, y_data, p, neg_r_inv, mu);

  // 模加减操作 - 优化分支
  u32 sum = gpu_stage3_mont_add(x, y, p);
  u32 diff = gpu_stage3_mont_sub(x, y, p);

  // 写回结果 - 合并写入
  a[idx1] = sum;
  a[idx2] = diff;
}

template <int MODE>
__global__ void ntt_stage3_kernel_inv(u32 *a, const u32 *twiddles_inv, int mid,
                                      int n, u32 p, u32 neg_r_inv, u64 mu) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total_butterflies = n >> 1;

  if (tid >= total_butterflies)
    return;

  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;

  // 预计算索引
  int tw_idx = (mid - 1) + k;
  int idx1 = j + k;
  int idx2 = idx1 + mid;

  // 寄存器重用
  u32 w_mont = twiddles_inv[tw_idx];
  u32 x = a[idx1];
  u32 y = a[idx2];

  // 逆蝶形操作
  u32 sum = gpu_stage3_mont_add(x, y, p);
  u32 temp = gpu_stage3_mont_sub(x, y, p);
  u32 diff = mod_mul_unified_stage3<MODE>(w_mont, temp, p, neg_r_inv, mu);

  a[idx1] = sum;
  a[idx2] = diff;
}

// 点乘核函数 - Stage 1优化
__global__ void pointwise_multiply_stage3(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                          u32 neg_r_inv) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    // 寄存器重用
    u32 a_val = a[idx];
    u32 b_val = b[idx];
    ab[idx] = gpu_stage3_mont_mul(a_val, b_val, p, neg_r_inv);
  }
}

// =============================================================================
// Stage 1优化主函数
// =============================================================================
template <int MODE>
void poly_multiply_ntt_gpu_stage3_core(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
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
  if constexpr (MODE == STAGE1_MUL_BARRETT) {
    mu = ((unsigned __int128)1 << 64) / p;
  }

  // 预计算旋转因子表
  std::vector<u32> tw_host(n_expanded);
  std::vector<u32> tw_inv_host(n_expanded);
  generate_twiddle_table_stage3(tw_host, n_expanded, p, omega, false);
  generate_twiddle_table_stage3(tw_inv_host, n_expanded, p, omega, true);

  // Stage 1优化：使用全局内存但优化访问模式
  u32 *tw_gpu, *tw_inv_gpu;
  CHECK_CUDA_STAGE3(cudaMalloc(&tw_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_STAGE3(cudaMalloc(&tw_inv_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_STAGE3(cudaMemcpy(tw_gpu, tw_host.data(), n_expanded * sizeof(u32),
                               cudaMemcpyHostToDevice));
  CHECK_CUDA_STAGE3(cudaMemcpy(tw_inv_gpu, tw_inv_host.data(),
                               n_expanded * sizeof(u32),
                               cudaMemcpyHostToDevice));

  // GPU内存分配
  u32 *a_gpu, *b_gpu, *ab_gpu;
  CHECK_CUDA_STAGE3(cudaMalloc(&a_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_STAGE3(cudaMalloc(&b_gpu, n_expanded * sizeof(u32)));
  CHECK_CUDA_STAGE3(cudaMalloc(&ab_gpu, n_expanded * sizeof(u32)));

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
  CHECK_CUDA_STAGE3(cudaMemcpy(a_gpu, a_mont_host, n_expanded * sizeof(u32),
                               cudaMemcpyHostToDevice));
  CHECK_CUDA_STAGE3(cudaMemcpy(b_gpu, b_mont_host, n_expanded * sizeof(u32),
                               cudaMemcpyHostToDevice));

  // Stage 3优化：动态选择最优线程块大小
  int BLOCK_SIZE = 256;
  if (n_expanded <= 16384) {
    BLOCK_SIZE = 128; // 小规模问题使用较小的block size
  } else if (n_expanded >= 131072) {
    BLOCK_SIZE = 512; // 大规模问题使用较大的block size
  }

  // 前向NTT变换
  for (u32 mid = 1; mid < n_expanded; mid <<= 1) {
    int total_butterflies = n_expanded >> 1;
    int grid_size = (total_butterflies + BLOCK_SIZE - 1) / BLOCK_SIZE;

    ntt_stage3_kernel_fwd<MODE><<<grid_size, BLOCK_SIZE>>>(
        a_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);
    ntt_stage3_kernel_fwd<MODE><<<grid_size, BLOCK_SIZE>>>(
        b_gpu, tw_gpu, mid, n_expanded, p, neg_r_inv, mu);

    CHECK_CUDA_STAGE3(cudaGetLastError());
  }
  CHECK_CUDA_STAGE3(cudaDeviceSynchronize());

  // 点乘操作
  {
    int grid_size = (n_expanded + BLOCK_SIZE - 1) / BLOCK_SIZE;
    pointwise_multiply_stage3<<<grid_size, BLOCK_SIZE>>>(
        a_gpu, b_gpu, ab_gpu, n_expanded, p, neg_r_inv);
    CHECK_CUDA_STAGE3(cudaGetLastError());
  }
  CHECK_CUDA_STAGE3(cudaDeviceSynchronize());

  // 逆NTT变换
  for (u32 mid = n_expanded >> 1; mid >= 1; mid >>= 1) {
    int total_butterflies = n_expanded >> 1;
    int grid_size = (total_butterflies + BLOCK_SIZE - 1) / BLOCK_SIZE;

    ntt_stage3_kernel_inv<MODE><<<grid_size, BLOCK_SIZE>>>(
        ab_gpu, tw_inv_gpu, mid, n_expanded, p, neg_r_inv, mu);

    CHECK_CUDA_STAGE3(cudaGetLastError());
  }
  CHECK_CUDA_STAGE3(cudaDeviceSynchronize());

  // 下载结果
  u32 *ab_mont_host = new u32[n_expanded];
  CHECK_CUDA_STAGE3(cudaMemcpy(ab_mont_host, ab_gpu, n_expanded * sizeof(u32),
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
  CHECK_CUDA_STAGE3(cudaFree(a_gpu));
  CHECK_CUDA_STAGE3(cudaFree(b_gpu));
  CHECK_CUDA_STAGE3(cudaFree(ab_gpu));
  CHECK_CUDA_STAGE3(cudaFree(tw_gpu));
  CHECK_CUDA_STAGE3(cudaFree(tw_inv_gpu));
}

// =============================================================================
// 公开接口 - Stage 1优化版本
// =============================================================================
void poly_multiply_ntt_gpu_stage3_mont(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                       u32 omega) {
  poly_multiply_ntt_gpu_stage3_core<STAGE1_MUL_MONT>(a, b, ab, n, p, omega);
}

void poly_multiply_ntt_gpu_stage3_naive(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                        u32 omega) {
  poly_multiply_ntt_gpu_stage3_core<STAGE1_MUL_NAIVE>(a, b, ab, n, p, omega);
}

void poly_multiply_ntt_gpu_stage3_barrett(u32 *a, u32 *b, u32 *ab, u32 n, u32 p,
                                          u32 omega) {
  poly_multiply_ntt_gpu_stage3_core<STAGE1_MUL_BARRETT>(a, b, ab, n, p, omega);
}