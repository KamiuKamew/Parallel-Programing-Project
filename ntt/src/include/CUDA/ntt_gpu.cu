#include "ntt.h"
#include <stdio.h>

// =============== 主机端模运算函数 ===============

u64 host_mod_mul(u64 a, u64 b, u64 mod) {
  return ((unsigned __int128)a * b) % mod;
}

u64 host_mod_pow(u64 base, u64 exp, u64 mod) {
  u64 result = 1;
  base %= mod;
  while (exp > 0) {
    if (exp & 1) {
      result = host_mod_mul(result, base, mod);
    }
    base = host_mod_mul(base, base, mod);
    exp >>= 1;
  }
  return result;
}

// =============== GPU设备函数实现 ===============

__device__ u64 gpu_mod_add(u64 a, u64 b, u64 mod) {
  u64 result = a + b;
  return (result >= mod) ? (result - mod) : result;
}

__device__ u64 gpu_mod_sub(u64 a, u64 b, u64 mod) {
  return (a >= b) ? (a - b) : (a + mod - b);
}

__device__ u64 gpu_mod_mul_naive(u64 a, u64 b, u64 mod) {
  // 基础模乘算法 - 使用128位中间结果避免溢出
  unsigned long long ah = a >> 32;
  unsigned long long al = a & 0xFFFFFFFF;
  unsigned long long bh = b >> 32;
  unsigned long long bl = b & 0xFFFFFFFF;

  unsigned long long p1 = al * bl;
  unsigned long long p2 = al * bh;
  unsigned long long p3 = ah * bl;
  unsigned long long p4 = ah * bh;

  unsigned long long carry =
      ((p1 >> 32) + (p2 & 0xFFFFFFFF) + (p3 & 0xFFFFFFFF)) >> 32;
  unsigned long long high = p4 + (p2 >> 32) + (p3 >> 32) + carry;
  unsigned long long low =
      p1 + ((p2 & 0xFFFFFFFF) << 32) + ((p3 & 0xFFFFFFFF) << 32);

  // 简化的模运算，这里使用朴素方法
  // 注意：对于大数取模，可能需要更精确的算法
  if (high == 0) {
    return low % mod;
  } else {
    // 对于高位不为0的情况，需要更复杂的处理
    // 这里先使用简化版本，后续可以优化
    return ((high % mod) * (0x100000000ULL % mod) + (low % mod)) % mod;
  }
}

__device__ u64 gpu_mod_pow(u64 base, u64 exp, u64 mod) {
  u64 result = 1;
  base %= mod;
  while (exp > 0) {
    if (exp & 1) {
      result = gpu_mod_mul_naive(result, base, mod);
    }
    base = gpu_mod_mul_naive(base, base, mod);
    exp >>= 1;
  }
  return result;
}

// =============== GPU Kernel函数实现 ===============

// NTT正变换kernel - 基于线程全局索引
__global__ void ntt_forward_gpu_kernel(u64 *a, u64 n, u64 mid, u64 p, u64 Wn) {
  // 计算全局线程索引
  u64 thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  u64 total_butterflies = n / 2; // 每个阶段的蝶形运算总数

  if (thread_id >= total_butterflies)
    return;

  // 计算该线程负责的蝶形运算的j和k
  u64 butterflies_per_group = mid; // 每组j循环的蝶形运算数
  u64 group_id = thread_id / butterflies_per_group;
  u64 k = thread_id % butterflies_per_group;
  u64 j = group_id * (mid << 1);

  // 计算旋转因子 w = Wn^k
  u64 w = gpu_mod_pow(Wn, k, p);

  // 蝶形运算
  u64 x = a[j + k];
  u64 y = gpu_mod_mul_naive(w, a[j + k + mid], p);

  a[j + k] = gpu_mod_add(x, y, p);
  a[j + k + mid] = gpu_mod_sub(x, y, p);
}

// NTT逆变换kernel - 基于线程全局索引
__global__ void ntt_inverse_gpu_kernel(u64 *a, u64 n, u64 mid, u64 p, u64 Wn) {
  // 计算全局线程索引
  u64 thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  u64 total_butterflies = n / 2; // 每个阶段的蝶形运算总数

  if (thread_id >= total_butterflies)
    return;

  // 计算该线程负责的蝶形运算的j和k
  u64 butterflies_per_group = mid; // 每组j循环的蝶形运算数
  u64 group_id = thread_id / butterflies_per_group;
  u64 k = thread_id % butterflies_per_group;
  u64 j = group_id * (mid << 1);

  // 计算旋转因子
  u64 w = gpu_mod_pow(Wn, k, p);

  // 蝶形运算（逆变换）
  u64 x = a[j + k];
  u64 y = a[j + k + mid];

  a[j + k] = gpu_mod_add(x, y, p);
  a[j + k + mid] = gpu_mod_sub(x, gpu_mod_mul_naive(w, y, p), p);
}

// =============== 主机端接口函数实现 ===============

void ntt_forward_gpu(u64 *a, u64 n, u64 p, u64 omega) {
  // 分配GPU内存
  u64 *d_a;
  CUDA_CHECK(cudaMalloc(&d_a, n * sizeof(u64)));
  CUDA_CHECK(cudaMemcpy(d_a, a, n * sizeof(u64), cudaMemcpyHostToDevice));

  // 执行NTT的每个阶段
  for (u64 mid = 1; mid < n; mid <<= 1) {
    u64 Wn = host_mod_pow(omega, (p - 1) / (mid << 1), p); // 在CPU上计算

    // 计算grid和block的大小，确保不超过CUDA限制
    u64 total_threads = n / 2; // 每个阶段处理n/2个蝶形运算
    u64 threads_per_block = (mid < 1024) ? mid : 1024; // 限制每个block的线程数
    u64 num_blocks =
        (total_threads + threads_per_block - 1) / threads_per_block;

    // 启动kernel
    ntt_forward_gpu_kernel<<<num_blocks, threads_per_block>>>(d_a, n, mid, p,
                                                              Wn);
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  // 将结果拷贝回主机
  CUDA_CHECK(cudaMemcpy(a, d_a, n * sizeof(u64), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaFree(d_a));
}

void ntt_inverse_gpu(u64 *a, u64 n, u64 p, u64 omega) {
  // 计算逆元
  u64 omega_inv = host_mod_pow(omega, p - 2, p); // 费马小定理求逆元

  // 直接调用正变换，但使用omega的逆元（与CPU版本一致）
  ntt_forward_gpu(a, n, p, omega_inv);

  // 乘以n的逆元
  u64 n_inv = host_mod_pow(n, p - 2, p);
  for (u64 i = 0; i < n; i++) {
    a[i] = host_mod_mul(a[i], n_inv, p);
  }
}

void poly_multiply_ntt_gpu(u64 *a, u64 *b, u64 *result, u64 n, u64 p,
                           u64 omega) {
  printf("DEBUG: GPU poly_multiply called with n=%lu\n", n);

  // 简化版本：仅用于验证GPU NTT核心功能
  // 扩展到2的幂
  u64 n_expanded = 1;
  while (n_expanded < 2 * n - 1) {
    n_expanded <<= 1;
  }

  printf("DEBUG: n_expanded=%lu\n", n_expanded);

  // 分配扩展数组并初始化
  u64 *a_expanded = new u64[n_expanded]();
  u64 *b_expanded = new u64[n_expanded]();

  // 复制原始数据，其余位置填充0
  for (u64 i = 0; i < n; i++) {
    a_expanded[i] = a[i];
    b_expanded[i] = b[i];
  }
  for (u64 i = n; i < n_expanded; i++) {
    a_expanded[i] = 0;
    b_expanded[i] = 0;
  }

  // GPU NTT前向变换
  ntt_forward_gpu(a_expanded, n_expanded, p, omega);
  ntt_forward_gpu(b_expanded, n_expanded, p, omega);

  // 点乘
  for (u64 i = 0; i < n_expanded; i++) {
    a_expanded[i] = host_mod_mul(a_expanded[i], b_expanded[i], p);
  }

  // GPU NTT逆变换
  ntt_inverse_gpu(a_expanded, n_expanded, p, omega);

  // 复制结果
  for (u64 i = 0; i < 2 * n - 1; i++) {
    result[i] = a_expanded[i];
  }

  delete[] a_expanded;
  delete[] b_expanded;
}