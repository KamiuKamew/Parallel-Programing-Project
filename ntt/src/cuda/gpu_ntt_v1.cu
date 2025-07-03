#include "gpu_ntt_common.h"
#include <cuda_runtime.h>
#include <iostream>

// 预计算旋转因子表
void precompute_twiddle_factors(u64 **d_twiddles, u64 n, u64 p, u64 omega) {
  std::cout << "[优化1] 预计算旋转因子表..." << std::endl;

  // 计算总的旋转因子数量：对于每个mid层级，需要mid个旋转因子
  u64 total_twiddles = 0;
  for (u64 mid = 1; mid < n; mid <<= 1) {
    total_twiddles += mid;
  }

  std::cout << "[优化1] 总共需要 " << total_twiddles << " 个旋转因子"
            << std::endl;

  // 在CPU上预计算所有旋转因子
  u64 *twiddles_cpu = new u64[total_twiddles];
  u64 offset = 0;

  for (u64 mid = 1; mid < n; mid <<= 1) {
    u64 w_len = mod_pow(omega, (p - 1) / (mid << 1), p);

    // 预计算这个mid层级的所有旋转因子
    for (u64 k = 0; k < mid; k++) {
      twiddles_cpu[offset + k] = mod_pow(w_len, k, p);
    }

    offset += mid;
  }

  // 将旋转因子复制到GPU
  cudaMalloc((void **)d_twiddles, total_twiddles * sizeof(u64));
  cudaMemcpy(*d_twiddles, twiddles_cpu, total_twiddles * sizeof(u64),
             cudaMemcpyHostToDevice);

  delete[] twiddles_cpu;
  std::cout << "[优化1] 旋转因子预计算完成" << std::endl;
}

// GPU kernel：并行位逆序排列
__global__ void bit_reverse_kernel_optimized(u64 *data, u64 n, int bits) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= n)
    return;

  int j = 0;
  for (int k = 0; k < bits; k++) {
    if (idx & (1 << k)) {
      j |= (1 << (bits - 1 - k));
    }
  }

  if (idx < j) {
    u64 temp = data[idx];
    data[idx] = data[j];
    data[j] = temp;
  }
}

// 优化版GPU kernel：使用预计算的旋转因子表
__global__ void butterfly_kernel_optimized(u64 *data, u64 n, u64 mid,
                                           u64 *twiddles, u64 twiddle_offset,
                                           u64 p) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int total_operations = n / 2;
  if (idx >= total_operations)
    return;

  // 计算(j,k)坐标
  int block_id = idx / mid;
  int k = idx % mid;
  int j = block_id * (mid << 1);

  int i1 = j + k;
  int i2 = j + k + mid;

  // 关键改进：直接从预计算表中获取旋转因子，消除O(mid)循环
  u64 w = twiddles[twiddle_offset + k];

  u64 u_val = data[i1];
  u64 v_val = mod_mul(data[i2], w, p);
  data[i1] = mod_add(u_val, v_val, p);
  data[i2] = mod_sub(u_val, v_val, p);
}

// 优化版GPU kernel：GPU归一化
__global__ void normalize_kernel(u64 *data, u64 n, u64 n_inv, u64 p) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    data[idx] = mod_mul(data[idx], n_inv, p);
  }
}

// 优化版GPU kernel：逐点乘法
__global__ void pointwise_mul_kernel_optimized(u64 *a, u64 *b, u64 *result,
                                               u64 n, u64 p) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    result[idx] = mod_mul(a[idx], b[idx], p);
  }
}

// 优化版GPU正变换
__host__ void ntt_forward_gpu_optimized(CudaU64Memory &d_data, u64 n, u64 p,
                                        u64 omega, u64 *d_twiddles) {
  std::cout << "[优化1] GPU正变换开始, n=" << n << std::endl;

  // 计算比特数
  int bits = 0;
  u64 temp = n - 1;
  while (temp) {
    bits++;
    temp >>= 1;
  }

  // 位逆序排列
  dim3 blockSize(256);
  dim3 gridSize((n + blockSize.x - 1) / blockSize.x);
  bit_reverse_kernel_optimized<<<gridSize, blockSize>>>(d_data.get(), n, bits);
  syncAndCheck("GPU位逆序排列");

  // 蝶形运算 - 使用预计算的旋转因子
  u64 twiddle_offset = 0;
  for (u64 mid = 1; mid < n; mid <<= 1) {
    int total_operations = n / 2;
    dim3 butterflyGrid((total_operations + blockSize.x - 1) / blockSize.x);

    butterfly_kernel_optimized<<<butterflyGrid, blockSize>>>(
        d_data.get(), n, mid, d_twiddles, twiddle_offset, p);
    syncAndCheck("GPU蝶形运算");

    twiddle_offset += mid;
  }

  std::cout << "[优化1] GPU正变换完成" << std::endl;
}

// 优化版GPU逆变换
__host__ void ntt_inverse_gpu_optimized(CudaU64Memory &d_data, u64 n, u64 p,
                                        u64 omega, u64 *d_twiddles_inv) {
  std::cout << "[优化1] GPU逆变换开始, n=" << n << std::endl;

  // 修正：逆变换 = 用正变换函数 + 逆变换旋转因子表
  // 关键修复：传递原始omega参数，使用预计算的逆变换旋转因子表
  ntt_forward_gpu_optimized(d_data, n, p, omega, d_twiddles_inv);

  // GPU归一化
  u64 n_inv = mod_inv(n, p);
  dim3 blockSize(256);
  dim3 gridSize((n + blockSize.x - 1) / blockSize.x);
  normalize_kernel<<<gridSize, blockSize>>>(d_data.get(), n, n_inv, p);
  syncAndCheck("GPU归一化");

  std::cout << "[优化1] GPU逆变换完成" << std::endl;
}

// 优化版多项式乘法
void poly_multiply_gpu_optimized_v1(u64 *a, u64 *b, u64 *result, u64 n, u64 p) {
  // 使用传入的模数参数，omega固定为3（题目保证所有模数的原根都是3）
  const u64 omega = 3;

  std::cout << "=== GPU优化版本1：预计算旋转因子 ===" << std::endl;
  std::cout << "参数: p=" << p << ", omega=" << omega << ", n=" << n
            << std::endl;

  // 计算扩展大小
  u64 n_expanded = 1;
  while (n_expanded < 2 * n - 1) {
    n_expanded <<= 1;
  }

  std::cout << "扩展规模: " << n << " -> " << n_expanded << std::endl;

  // 验证参数有效性
  if ((p - 1) % n_expanded != 0) {
    std::cout << "❌ 错误: (p-1)=" << (p - 1)
              << "不能被n_expanded=" << n_expanded << "整除" << std::endl;
    return;
  }

  try {
    // CPU内存分配
    CpuU64Memory a_exp(n_expanded);
    CpuU64Memory b_exp(n_expanded);
    CpuU64Memory result_exp(n_expanded);

    // 扩展数组
    for (u64 i = 0; i < n; i++) {
      a_exp[i] = a[i];
      b_exp[i] = b[i];
    }
    for (u64 i = n; i < n_expanded; i++) {
      a_exp[i] = 0;
      b_exp[i] = 0;
    }

    // GPU内存分配
    CudaU64Memory d_a(n_expanded);
    CudaU64Memory d_b(n_expanded);
    CudaU64Memory d_result(n_expanded);

    // 数据传输到GPU
    d_a.copyFromHost(a_exp.get(), n_expanded);
    d_b.copyFromHost(b_exp.get(), n_expanded);

    // 预计算旋转因子
    u64 *d_twiddles;
    u64 *d_twiddles_inv;
    precompute_twiddle_factors(&d_twiddles, n_expanded, p, omega);

    u64 omega_inv = mod_inv(omega, p);
    precompute_twiddle_factors(&d_twiddles_inv, n_expanded, p, omega_inv);

    // 正变换
    ntt_forward_gpu_optimized(d_a, n_expanded, p, omega, d_twiddles);
    ntt_forward_gpu_optimized(d_b, n_expanded, p, omega, d_twiddles);

    // 逐点乘法
    dim3 blockSize(256);
    dim3 gridSize((n_expanded + blockSize.x - 1) / blockSize.x);
    pointwise_mul_kernel_optimized<<<gridSize, blockSize>>>(
        d_a.get(), d_b.get(), d_result.get(), n_expanded, p);
    syncAndCheck("GPU逐点乘法");

    // 逆变换
    ntt_inverse_gpu_optimized(d_result, n_expanded, p, omega, d_twiddles_inv);

    // 结果传输回CPU
    d_result.copyToHost(result_exp.get(), n_expanded);

    // 复制结果
    for (u64 i = 0; i < 2 * n - 1; i++) {
      result[i] = result_exp[i];
    }

    // 清理GPU内存
    cudaFree(d_twiddles);
    cudaFree(d_twiddles_inv);

    std::cout << "[优化1] 多项式乘法完成！" << std::endl;

  } catch (const std::exception &e) {
    std::cout << "❌ 异常: " << e.what() << std::endl;
  }
}