#include "ntt.h"
#include <cassert>
#include <iostream>

// CUDA错误检查辅助函数
void check_cuda_error_serial(const char *msg) {
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::cerr << "CUDA Error (" << msg << "): " << cudaGetErrorString(err)
              << std::endl;
    exit(EXIT_FAILURE);
  }
}

// 在GPU设备上实现Montgomery模运算（完整版本）
__device__ u64 mont_mul_device(u64 a, u64 b, u64 p, u64 neg_r_inv) {
  // Montgomery乘法：计算 a * b * R^(-1) mod p
  // 其中 R = 2^64
  unsigned __int128 t = (unsigned __int128)a * b;
  u64 m = (u64)t * neg_r_inv; // (t mod R) * neg_r_inv mod R
  unsigned __int128 tmp = t + (unsigned __int128)m * p;
  u64 res = (u64)(tmp >> 64); // tmp / R
  return res >= p ? res - p : res;
}

__device__ u64 mont_add_device(u64 a, u64 b, u64 p) {
  return (a + b >= p) ? (a + b - p) : (a + b);
}

__device__ u64 mont_sub_device(u64 a, u64 b, u64 p) {
  return (a >= b) ? (a - b) : (a + p - b);
}

// 在GPU上计算Montgomery域的快速幂
__device__ u64 mont_pow_device(u64 base, u64 exp, u64 p, u64 neg_r_inv,
                               u64 r2) {
  u64 result = mont_mul_device(1ULL, r2, p, neg_r_inv); // Montgomery形式的1
  while (exp > 0) {
    if (exp & 1) {
      result = mont_mul_device(result, base, p, neg_r_inv);
    }
    base = mont_mul_device(base, base, p, neg_r_inv);
    exp >>= 1;
  }
  return result;
}

// 将普通数转换为Montgomery形式
__device__ u64 to_mont_device(u64 a, u64 p, u64 r2, u64 neg_r_inv) {
  return mont_mul_device(a, r2, p, neg_r_inv);
}

// 将Montgomery形式转换为普通数
__device__ u64 from_mont_device(u64 a_mont, u64 p, u64 neg_r_inv) {
  return mont_mul_device(a_mont, 1ULL, p, neg_r_inv);
}

// 第二步改进：预处理所有旋转因子
// 为每个(mid, k)组合预计算 w_k = omega^((p-1)/(2*mid) * k)
void precompute_twiddle_factors(u64 **d_twiddles, u64 n, u64 p, u64 omega_mont,
                                u64 neg_r_inv, u64 r2) {
  std::cout << "[第二步改进] 预处理旋转因子..." << std::endl;

  // 计算总的旋转因子数量：对于每个mid层级，需要mid个旋转因子
  u64 total_twiddles = 0;
  for (u64 mid = 1; mid < n; mid <<= 1) {
    total_twiddles += mid;
  }

  // 在CPU上计算所有旋转因子
  u64 *twiddles_cpu = new u64[total_twiddles];
  u64 offset = 0;

  MontMod<u64> mont_mod(p);

  for (u64 mid = 1; mid < n; mid <<= 1) {
    // 计算 Wn = omega^((p-1)/(2*mid))
    u64 exp = (p - 1) / (mid << 1);
    u64 Wn_mont =
        mont_mod.pow(mont_mod.from_T(3), exp); // 使用CPU Montgomery计算

    // 为这个mid层级预计算所有旋转因子
    u64 w_mont = mont_mod.from_T(1); // Montgomery形式的1
    for (u64 k = 0; k < mid; ++k) {
      twiddles_cpu[offset + k] = w_mont;
      w_mont = mont_mod.mul(w_mont, Wn_mont); // w *= Wn
    }
    offset += mid;
  }

  // 将旋转因子复制到GPU
  cudaMalloc((void **)d_twiddles, total_twiddles * sizeof(u64));
  cudaMemcpy(*d_twiddles, twiddles_cpu, total_twiddles * sizeof(u64),
             cudaMemcpyHostToDevice);
  check_cuda_error_serial("precompute twiddle factors");

  delete[] twiddles_cpu;
  std::cout << "[第二步改进] 预处理了 " << total_twiddles << " 个旋转因子"
            << std::endl;
}

// 第二步改进：重构的NTT正变换kernel（便于并行化但仍串行执行）
__global__ void ntt_forward_parallel_ready_kernel(u64 *a_mont, u64 n, u64 p,
                                                  u64 *d_twiddles,
                                                  u64 neg_r_inv) {
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    // 第二步重构：使用预处理的旋转因子，重构循环结构
    u64 twiddle_offset = 0; // 当前层级在旋转因子数组中的偏移

    for (u64 mid = 1; mid < n; mid <<= 1) {
      // 重构的关键：将(j,k)组合重新组织便于并行
      // 在第二步中仍用串行执行，但结构便于第三步并行化

      u64 pairs_count = n / 2; // 每个mid层级的蝶形运算对数

      // 重构循环：将原来的嵌套循环展开为单层循环
      for (u64 tid = 0; tid < pairs_count; ++tid) {
        // 从线程ID反推出(j,k)坐标
        u64 k = tid % mid;            // k在[0, mid)范围内
        u64 j_group = tid / mid;      // j的组索引
        u64 j = j_group * (mid << 1); // 实际的j值

        // 使用预处理的旋转因子（无需重复计算pow）
        u64 w_mont = d_twiddles[twiddle_offset + k];

        // 蝶形运算（逻辑完全一致）
        u64 x_mont = a_mont[j + k];
        u64 y_mont = mont_mul_device(w_mont, a_mont[j + k + mid], p, neg_r_inv);

        a_mont[j + k] = mont_add_device(x_mont, y_mont, p);
        a_mont[j + k + mid] = mont_sub_device(x_mont, y_mont, p);
      }

      twiddle_offset += mid; // 移动到下一层级的旋转因子
    }
  }
}

// 第二步改进：重构的NTT逆变换kernel（便于并行化但仍串行执行）
__global__ void ntt_inverse_parallel_ready_kernel(u64 *a_mont, u64 n, u64 p,
                                                  u64 *d_twiddles_inv,
                                                  u64 neg_r_inv, u64 r2) {
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    // 第二步重构：使用预处理的逆变换旋转因子
    u64 twiddle_offset = 0;

    for (u64 mid = n >> 1; mid > 0; mid >>= 1) {
      u64 pairs_count = n / 2;

      // 重构循环：将原来的嵌套循环展开为单层循环
      for (u64 tid = 0; tid < pairs_count; ++tid) {
        u64 k = tid % mid;
        u64 j_group = tid / mid;
        u64 j = j_group * (mid << 1);

        // 使用预处理的逆变换旋转因子
        u64 w_mont = d_twiddles_inv[twiddle_offset + k];

        // 逆蝶形运算（逻辑完全一致）
        u64 x_mont = a_mont[j + k];
        u64 y_mont = a_mont[j + k + mid];

        a_mont[j + k] = mont_add_device(x_mont, y_mont, p);
        a_mont[j + k + mid] = mont_mul_device(
            w_mont, mont_sub_device(x_mont, y_mont, p), p, neg_r_inv);
      }

      twiddle_offset += mid;
    }

    // 乘以n的逆元
    u64 n_mont = to_mont_device(n, p, r2, neg_r_inv);
    u64 inv_n_mont = mont_pow_device(n_mont, p - 2, p, neg_r_inv, r2);
    for (u64 i = 0; i < n; ++i) {
      a_mont[i] = mont_mul_device(a_mont[i], inv_n_mont, p, neg_r_inv);
    }
  }
}

// Montgomery域的逐点乘法
__global__ void pointwise_mul_serial_kernel(u64 *a_mont, u64 *b_mont,
                                            u64 *ab_mont, u64 n, u64 p,
                                            u64 neg_r_inv) {
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    for (u64 i = 0; i < n; ++i) {
      ab_mont[i] = mont_mul_device(a_mont[i], b_mont[i], p, neg_r_inv);
    }
  }
}

// 第一步：真正的基础串行版本 - Montgomery NTT算法
void poly_multiply_ntt_cuda_serial(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                   u64 OMEGA) {
  std::cout << "[CUDA Serial Montgomery] 开始真正的串行Montgomery NTT..."
            << std::endl;

  // 计算扩展后的长度
  u64 n_expanded = expand_n(2 * n - 1);

  // 在CPU上扩展数组
  u64 *a_expanded = expand_a(a, n, n_expanded);
  u64 *b_expanded = expand_a(b, n, n_expanded);

  // 位逆序排列
  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  // 计算Montgomery参数（复制CPU上MontMod的逻辑）
  // 计算 R^2 mod p，其中 R = 2^64
  unsigned __int128 r_val_mod_n = 1;
  for (int i = 0; i < 64; ++i) {
    r_val_mod_n = (r_val_mod_n << 1) % p;
  }
  u64 r2 = (unsigned __int128)r_val_mod_n * r_val_mod_n % p;

  // 计算 -p^(-1) mod R 使用牛顿迭代法
  u64 inv = 1;
  for (int i = 0; i < 6; ++i) { // 6次迭代对于64位足够
    inv = (u64)((unsigned __int128)inv * (2 - (unsigned __int128)p * inv));
  }
  u64 neg_r_inv = -inv;

  // 使用CPU上的Montgomery乘法转换omega到Montgomery域
  MontMod<u64> mont_mod(p);
  u64 omega_mont = mont_mod.from_T(OMEGA);

  // 计算omega的逆元并转换到Montgomery域
  Mod<u64> mod(p);
  u64 omega_inv = mod.inv(OMEGA);
  u64 omega_inv_mont = mont_mod.from_T(omega_inv);

  // 分配GPU内存
  u64 *d_a_mont, *d_b_mont, *d_ab_mont;
  cudaMalloc((void **)&d_a_mont, n_expanded * sizeof(u64));
  cudaMalloc((void **)&d_b_mont, n_expanded * sizeof(u64));
  cudaMalloc((void **)&d_ab_mont, n_expanded * sizeof(u64));
  check_cuda_error_serial("GPU memory allocation");

  // 转换数据到Montgomery域并复制到GPU
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
  check_cuda_error_serial("copy to GPU");

  // 在GPU上执行串行NTT正变换
  ntt_forward_serial_kernel<<<1, 1>>>(d_a_mont, n_expanded, p, omega_mont,
                                      neg_r_inv, r2);
  cudaDeviceSynchronize();
  check_cuda_error_serial("ntt_forward_serial_kernel");

  ntt_forward_serial_kernel<<<1, 1>>>(d_b_mont, n_expanded, p, omega_mont,
                                      neg_r_inv, r2);
  cudaDeviceSynchronize();
  check_cuda_error_serial("ntt_forward_serial_kernel");

  // 逐点乘法
  pointwise_mul_serial_kernel<<<1, 1>>>(d_a_mont, d_b_mont, d_ab_mont,
                                        n_expanded, p, neg_r_inv);
  cudaDeviceSynchronize();
  check_cuda_error_serial("pointwise_mul_serial_kernel");

  // 逆变换
  ntt_inverse_serial_kernel<<<1, 1>>>(d_ab_mont, n_expanded, p, omega_inv_mont,
                                      neg_r_inv, r2);
  cudaDeviceSynchronize();
  check_cuda_error_serial("ntt_inverse_serial_kernel");

  // 将结果复制回CPU并转换出Montgomery域
  u64 *result_mont_cpu = new u64[n_expanded];
  cudaMemcpy(result_mont_cpu, d_ab_mont, n_expanded * sizeof(u64),
             cudaMemcpyDeviceToHost);
  check_cuda_error_serial("copy to CPU");

  u64 *result_expanded = new u64[n_expanded];
  for (u64 i = 0; i < n_expanded; ++i) {
    result_expanded[i] = mont_mod.to_T(result_mont_cpu[i]);
  }

  // 位逆序排列
  bit_reverse_permute(result_expanded, n_expanded);

  // 复制结果
  for (u64 i = 0; i < 2 * n - 1; ++i) {
    ab[i] = result_expanded[i];
  }

  // 清理内存
  cudaFree(d_a_mont);
  cudaFree(d_b_mont);
  cudaFree(d_ab_mont);
  delete[] a_expanded;
  delete[] b_expanded;
  delete[] a_mont_cpu;
  delete[] b_mont_cpu;
  delete[] result_mont_cpu;
  delete[] result_expanded;

  std::cout << "[CUDA Serial Montgomery] 真正的串行Montgomery NTT完成"
            << std::endl;
}

// 第二步：便于并行化的串行版本 - 在第一步基础上改进
void poly_multiply_ntt_cuda_parallel_ready(u64 *a, u64 *b, u64 *ab, u64 n,
                                           u64 p, u64 OMEGA) {
  std::cout << "[CUDA 第二步] 便于并行化的串行版本..." << std::endl;

  // 基础设置（与第一步相同）
  u64 n_expanded = expand_n(2 * n - 1);
  u64 *a_expanded = expand_a(a, n, n_expanded);
  u64 *b_expanded = expand_a(b, n, n_expanded);
  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  // Montgomery参数计算（与第一步相同）
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

  // 第二步改进：预处理旋转因子
  u64 *d_twiddles, *d_twiddles_inv;
  precompute_twiddle_factors(&d_twiddles, n_expanded, p, omega_mont, neg_r_inv,
                             r2);
  precompute_twiddle_factors(&d_twiddles_inv, n_expanded, p, omega_inv_mont,
                             neg_r_inv, r2);

  // GPU内存分配和数据传输（与第一步相同）
  u64 *d_a_mont, *d_b_mont, *d_ab_mont;
  cudaMalloc((void **)&d_a_mont, n_expanded * sizeof(u64));
  cudaMalloc((void **)&d_b_mont, n_expanded * sizeof(u64));
  cudaMalloc((void **)&d_ab_mont, n_expanded * sizeof(u64));
  check_cuda_error_serial("GPU memory allocation");

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
  check_cuda_error_serial("copy to GPU");

  // 第二步改进：使用重构的kernel
  ntt_forward_parallel_ready_kernel<<<1, 1>>>(d_a_mont, n_expanded, p,
                                              d_twiddles, neg_r_inv);
  cudaDeviceSynchronize();
  check_cuda_error_serial("ntt_forward_parallel_ready_kernel");

  ntt_forward_parallel_ready_kernel<<<1, 1>>>(d_b_mont, n_expanded, p,
                                              d_twiddles, neg_r_inv);
  cudaDeviceSynchronize();
  check_cuda_error_serial("ntt_forward_parallel_ready_kernel");

  // 逐点乘法（与第一步相同）
  pointwise_mul_serial_kernel<<<1, 1>>>(d_a_mont, d_b_mont, d_ab_mont,
                                        n_expanded, p, neg_r_inv);
  cudaDeviceSynchronize();
  check_cuda_error_serial("pointwise_mul_serial_kernel");

  // 第二步改进：使用重构的逆变换kernel
  ntt_inverse_parallel_ready_kernel<<<1, 1>>>(d_ab_mont, n_expanded, p,
                                              d_twiddles_inv, neg_r_inv, r2);
  cudaDeviceSynchronize();
  check_cuda_error_serial("ntt_inverse_parallel_ready_kernel");

  // 结果处理（与第一步相同）
  u64 *result_mont_cpu = new u64[n_expanded];
  cudaMemcpy(result_mont_cpu, d_ab_mont, n_expanded * sizeof(u64),
             cudaMemcpyDeviceToHost);
  check_cuda_error_serial("copy to CPU");

  u64 *result_expanded = new u64[n_expanded];
  for (u64 i = 0; i < n_expanded; ++i) {
    result_expanded[i] = mont_mod.to_T(result_mont_cpu[i]);
  }
  bit_reverse_permute(result_expanded, n_expanded);
  for (u64 i = 0; i < 2 * n - 1; ++i) {
    ab[i] = result_expanded[i];
  }

  // 清理内存（包括新的旋转因子内存）
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

  std::cout << "[CUDA 第二步] 便于并行化的串行版本完成" << std::endl;
}