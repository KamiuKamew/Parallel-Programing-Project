#include "ntt.h"
#include <cassert>
#include <iostream>

// CUDA错误检查辅助函数
void check_cuda_error_step3(const char *msg) {
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::cerr << "CUDA Error (" << msg << "): " << cudaGetErrorString(err)
              << std::endl;
    exit(EXIT_FAILURE);
  }
}

// 在GPU设备上实现Montgomery模运算（与前面步骤相同）
__device__ u64 mont_mul_device(u64 a, u64 b, u64 p, u64 neg_r_inv) {
  unsigned __int128 t = (unsigned __int128)a * b;
  u64 m = (u64)t * neg_r_inv;
  unsigned __int128 tmp = t + (unsigned __int128)m * p;
  u64 res = (u64)(tmp >> 64);
  return res >= p ? res - p : res;
}

__device__ u64 mont_add_device(u64 a, u64 b, u64 p) {
  return (a + b >= p) ? (a + b - p) : (a + b);
}

__device__ u64 mont_sub_device(u64 a, u64 b, u64 p) {
  return (a >= b) ? (a - b) : (a + p - b);
}

__device__ u64 mont_pow_device(u64 base, u64 exp, u64 p, u64 neg_r_inv,
                               u64 r2) {
  u64 result = mont_mul_device(1ULL, r2, p, neg_r_inv);
  while (exp > 0) {
    if (exp & 1) {
      result = mont_mul_device(result, base, p, neg_r_inv);
    }
    base = mont_mul_device(base, base, p, neg_r_inv);
    exp >>= 1;
  }
  return result;
}

__device__ u64 to_mont_device(u64 a, u64 p, u64 r2, u64 neg_r_inv) {
  return mont_mul_device(a, r2, p, neg_r_inv);
}

// 第三步改进：预处理旋转因子（与第二步相同）
void precompute_twiddle_factors_step3(u64 **d_twiddles, u64 n, u64 p,
                                      u64 omega_mont) {
  std::cout << "[第三步并行] 预处理旋转因子..." << std::endl;

  u64 total_twiddles = 0;
  for (u64 mid = 1; mid < n; mid <<= 1) {
    total_twiddles += mid;
  }

  u64 *twiddles_cpu = new u64[total_twiddles];
  u64 offset = 0;

  MontMod<u64> mont_mod(p);

  for (u64 mid = 1; mid < n; mid <<= 1) {
    u64 exp = (p - 1) / (mid << 1);
    u64 Wn_mont = mont_mod.pow(omega_mont, exp);

    u64 w_mont = mont_mod.from_T(1);
    for (u64 k = 0; k < mid; ++k) {
      twiddles_cpu[offset + k] = w_mont;
      w_mont = mont_mod.mul(w_mont, Wn_mont);
    }
    offset += mid;
  }

  cudaMalloc((void **)d_twiddles, total_twiddles * sizeof(u64));
  cudaMemcpy(*d_twiddles, twiddles_cpu, total_twiddles * sizeof(u64),
             cudaMemcpyHostToDevice);
  check_cuda_error_step3("precompute twiddle factors");

  delete[] twiddles_cpu;
  std::cout << "[第三步并行] 预处理了 " << total_twiddles << " 个旋转因子"
            << std::endl;
}

// 第三步：真正的GPU并行NTT正变换kernel
__global__ void ntt_forward_parallel_kernel(u64 *a_mont, u64 n, u64 p,
                                            u64 *d_twiddles, u64 neg_r_inv,
                                            u64 mid, u64 twiddle_offset) {
  // 关键改进：将串行循环替换为GPU线程并行
  u64 tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < n / 2) {
    // 完全相同的线程ID到坐标映射逻辑（来自第二步）
    u64 block_id = tid / mid;
    u64 k = tid % mid;
    u64 j = block_id * (mid << 1);

    u64 w_mont = d_twiddles[twiddle_offset + k];

    u64 x_mont = a_mont[j + k];
    u64 y_mont = mont_mul_device(w_mont, a_mont[j + k + mid], p, neg_r_inv);

    a_mont[j + k] = mont_add_device(x_mont, y_mont, p);
    a_mont[j + k + mid] = mont_sub_device(x_mont, y_mont, p);
  }
}

// 第三步：真正的GPU并行NTT逆变换kernel
__global__ void ntt_inverse_parallel_kernel(u64 *a_mont, u64 n, u64 p,
                                            u64 *d_twiddles_inv, u64 neg_r_inv,
                                            u64 mid, u64 twiddle_offset) {
  // 关键改进：将串行循环替换为GPU线程并行
  u64 tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < n / 2) {
    // 完全相同的线程ID到坐标映射逻辑（来自第二步）
    u64 block_id = tid / mid;
    u64 k = tid % mid;
    u64 j = block_id * (mid << 1);

    u64 w_mont = d_twiddles_inv[twiddle_offset + k];

    u64 x_mont = a_mont[j + k];
    u64 y_mont = a_mont[j + k + mid];

    a_mont[j + k] = mont_add_device(x_mont, y_mont, p);
    a_mont[j + k + mid] = mont_mul_device(
        w_mont, mont_sub_device(x_mont, y_mont, p), p, neg_r_inv);
  }
}

// 第三步：真正的GPU并行逐点乘法kernel
__global__ void pointwise_mul_parallel_kernel(u64 *a_mont, u64 *b_mont,
                                              u64 *ab_mont, u64 n, u64 p,
                                              u64 neg_r_inv) {
  // 关键改进：将串行循环替换为GPU线程并行
  u64 i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < n) {
    ab_mont[i] = mont_mul_device(a_mont[i], b_mont[i], p, neg_r_inv);
  }
}

// 第三步：GPU并行归一化kernel（乘以逆元）
__global__ void normalize_parallel_kernel(u64 *a_mont, u64 n, u64 inv_n_mont,
                                          u64 p, u64 neg_r_inv) {
  u64 i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    a_mont[i] = mont_mul_device(a_mont[i], inv_n_mont, p, neg_r_inv);
  }
}

// 第三步：真正的GPU并行版本
void poly_multiply_ntt_cuda_parallel(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                     u64 OMEGA) {
  std::cout << "[CUDA 第三步] 真正的GPU并行版本..." << std::endl;

  // 基础设置（与前面步骤相同）
  u64 n_expanded = expand_n(2 * n - 1);
  u64 *a_expanded = expand_a(a, n, n_expanded);
  u64 *b_expanded = expand_a(b, n, n_expanded);
  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  // Montgomery参数计算（与前面步骤相同）
  unsigned __int128 r_val_mod_n = 1;
  for (int i = 0; i < 64; ++i) {
    r_val_mod_n = (r_val_mod_n << 1) % p;
  }
  u64 r2 = (unsigned __int128)r_val_mod_n * r_val_mod_n % p;
  u64 inv = 1;
  for (int i = 0; i < 6; ++i) {
    inv = (u64)((unsigned __int128)inv * (2 - (unsigned __int128)p * inv));
  }
  u64 neg_r_inv = -inv;

  MontMod<u64> mont_mod(p);
  u64 omega_mont = mont_mod.from_T(OMEGA);
  Mod<u64> mod(p);
  u64 omega_inv = mod.inv(OMEGA);
  u64 omega_inv_mont = mont_mod.from_T(omega_inv);

  // 第三步改进：预处理旋转因子
  u64 *d_twiddles, *d_twiddles_inv;
  precompute_twiddle_factors_step3(&d_twiddles, n_expanded, p, omega_mont);
  precompute_twiddle_factors_step3(&d_twiddles_inv, n_expanded, p,
                                   omega_inv_mont);

  // GPU内存分配和数据传输（与前面步骤相同）
  u64 *d_a_mont, *d_b_mont, *d_ab_mont;
  cudaMalloc((void **)&d_a_mont, n_expanded * sizeof(u64));
  cudaMalloc((void **)&d_b_mont, n_expanded * sizeof(u64));
  cudaMalloc((void **)&d_ab_mont, n_expanded * sizeof(u64));
  check_cuda_error_step3("GPU memory allocation");

  u64 *a_mont_cpu = new u64[n_expanded];
  u64 *b_mont_cpu = new u64[n_expanded];
  for (u64 i = 0; i < n_expanded; ++i) {
    a_mont_cpu[i] = mont_mod.from_T(a_expanded[i]);
    b_mont_cpu[i] = mont_mod.from_T(b_expanded[i]);
  }

  cudaMemcpy(d_a_mont, a_mont_cpu, n_expanded * sizeof(u64),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_b_mont, b_mont_cpu, n_expanded * sizeof(u64),
             cudaMemcpyHostToDevice);
  check_cuda_error_step3("copy to GPU");

  // 关键改进：GPU并行配置
  const int threads_per_block = 256;
  int blocks_for_ntt =
      (n_expanded / 2 + threads_per_block - 1) / threads_per_block;
  int blocks_for_pointwise =
      (n_expanded + threads_per_block - 1) / threads_per_block;

  std::cout << "[第三步并行] GPU配置: " << blocks_for_ntt << " blocks × "
            << threads_per_block
            << " threads = " << blocks_for_ntt * threads_per_block
            << " 总线程 (NTT需要 " << n_expanded / 2 << " 线程)" << std::endl;

  // 第三步改进：并行NTT正变换
  u64 twiddle_offset = 0;
  for (u64 mid = 1; mid < n_expanded; mid <<= 1) {
    ntt_forward_parallel_kernel<<<blocks_for_ntt, threads_per_block>>>(
        d_a_mont, n_expanded, p, d_twiddles, neg_r_inv, mid, twiddle_offset);
    cudaDeviceSynchronize();
    check_cuda_error_step3("ntt_forward_parallel_kernel");

    ntt_forward_parallel_kernel<<<blocks_for_ntt, threads_per_block>>>(
        d_b_mont, n_expanded, p, d_twiddles, neg_r_inv, mid, twiddle_offset);
    cudaDeviceSynchronize();
    check_cuda_error_step3("ntt_forward_parallel_kernel");

    twiddle_offset += mid;
  }

  // 第三步改进：并行逐点乘法
  pointwise_mul_parallel_kernel<<<blocks_for_pointwise, threads_per_block>>>(
      d_a_mont, d_b_mont, d_ab_mont, n_expanded, p, neg_r_inv);
  cudaDeviceSynchronize();
  check_cuda_error_step3("pointwise_mul_parallel_kernel");

  // 第三步改进：并行NTT逆变换
  for (u64 mid = n_expanded >> 1; mid > 0; mid >>= 1) {
    // 计算当前mid对应的正确旋转因子偏移量
    u64 twiddle_offset = 0;
    for (u64 stored_mid = 1; stored_mid < mid; stored_mid <<= 1) {
      twiddle_offset += stored_mid;
    }

    ntt_inverse_parallel_kernel<<<blocks_for_ntt, threads_per_block>>>(
        d_ab_mont, n_expanded, p, d_twiddles_inv, neg_r_inv, mid,
        twiddle_offset);
    cudaDeviceSynchronize();
    check_cuda_error_step3("ntt_inverse_parallel_kernel");
  }

  // 乘以n的逆元（使用并行GPU归一化）
  u64 n_mont = mont_mod.from_T(n_expanded);
  u64 inv_n_mont = mont_mod.inv(n_mont);

  // 使用专门的归一化kernel
  normalize_parallel_kernel<<<blocks_for_pointwise, threads_per_block>>>(
      d_ab_mont, n_expanded, inv_n_mont, p, neg_r_inv);
  cudaDeviceSynchronize();
  check_cuda_error_step3("normalize_parallel_kernel");

  // 结果处理（与前面步骤相同）
  u64 *result_mont_cpu = new u64[n_expanded];
  cudaMemcpy(result_mont_cpu, d_ab_mont, n_expanded * sizeof(u64),
             cudaMemcpyDeviceToHost);
  check_cuda_error_step3("copy to CPU");

  u64 *result_expanded = new u64[n_expanded];
  for (u64 i = 0; i < n_expanded; ++i) {
    result_expanded[i] = mont_mod.to_T(result_mont_cpu[i]);
  }
  bit_reverse_permute(result_expanded, n_expanded);
  for (u64 i = 0; i < 2 * n - 1; ++i) {
    ab[i] = result_expanded[i];
  }

  // 清理内存
  cudaFree(d_a_mont);
  cudaFree(d_b_mont);
  cudaFree(d_ab_mont);
  cudaFree(d_twiddles);
  cudaFree(d_twiddles_inv);
  delete[] a_expanded;
  delete[] b_expanded;
  delete[] a_mont_cpu;
  delete[] b_mont_cpu;
  delete[] result_mont_cpu;
  delete[] result_expanded;

  std::cout << "[CUDA 第三步] 真正的GPU并行版本完成" << std::endl;
}