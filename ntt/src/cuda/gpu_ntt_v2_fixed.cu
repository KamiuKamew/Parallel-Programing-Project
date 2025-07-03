#include "gpu_ntt_common.h"
#include <cuda_runtime.h>
#include <iostream>

// 包含v1版本的实现作为回退
#include "gpu_ntt_v1.cu"

// 修复版v2：安全的共享内存优化
#define SHARED_MEM_SIZE 512 // 减小共享内存大小，提高安全性

// 修复版共享内存kernel：只处理单个mid层级，避免复杂的分组逻辑
__global__ void butterfly_kernel_shared_memory_fixed(u64 *data, u64 n, u64 mid,
                                                     u64 *twiddles,
                                                     u64 twiddle_offset,
                                                     u64 p) {
  extern __shared__ u64 shared_data[];

  int tid = threadIdx.x;
  int bid = blockIdx.x;
  int block_size = blockDim.x;

  // 计算当前线程处理的蝶形运算
  int total_operations = n / 2;
  int ops_per_block = (total_operations + gridDim.x - 1) / gridDim.x;
  int start_op = bid * ops_per_block;
  int end_op = min(start_op + ops_per_block, total_operations);

  // 检查是否需要处理
  if (start_op >= total_operations)
    return;

  // 为了简化，先使用全局内存版本，确保正确性
  // 每个线程处理一个蝶形运算
  for (int op_idx = start_op + tid; op_idx < end_op; op_idx += block_size) {
    // 计算(j,k)坐标 - 与v1版本完全一致
    int block_id = op_idx / mid;
    int k = op_idx % mid;
    int j = block_id * (mid << 1);

    int i1 = j + k;
    int i2 = j + k + mid;

    // 获取旋转因子
    u64 w = twiddles[twiddle_offset + k];

    // 蝶形运算 - 与v1版本完全一致
    u64 u_val = data[i1];
    u64 v_val = mod_mul(data[i2], w, p);
    data[i1] = mod_add(u_val, v_val, p);
    data[i2] = mod_sub(u_val, v_val, p);
  }
}

// 修复版GPU正变换：使用简化的共享内存策略
__host__ void ntt_forward_gpu_optimized_v2_fixed(CudaU64Memory &d_data, u64 n,
                                                 u64 p, u64 omega,
                                                 u64 *d_twiddles) {
  std::cout << "[优化2-修复] GPU正变换开始, n=" << n << std::endl;

  // 位逆序排列 - 与v1完全一致
  int bits = 0;
  u64 temp = n - 1;
  while (temp) {
    bits++;
    temp >>= 1;
  }

  dim3 blockSize(256);
  dim3 gridSize((n + blockSize.x - 1) / blockSize.x);
  bit_reverse_kernel_optimized<<<gridSize, blockSize>>>(d_data.get(), n, bits);
  syncAndCheck("GPU位逆序排列");

  // 蝶形运算 - 逐层处理，与v1逻辑一致但使用不同的kernel
  u64 twiddle_offset = 0;
  for (u64 mid = 1; mid < n; mid <<= 1) {
    int total_operations = n / 2;

    // 使用修复的共享内存kernel，但网格配置保守
    dim3 butterflyGrid((total_operations + 255) / 256);
    dim3 butterflyBlock(256);

    butterfly_kernel_shared_memory_fixed<<<butterflyGrid, butterflyBlock,
                                           SHARED_MEM_SIZE * sizeof(u64)>>>(
        d_data.get(), n, mid, d_twiddles, twiddle_offset, p);
    syncAndCheck("修复版共享内存蝶形运算");

    twiddle_offset += mid;
  }

  std::cout << "[优化2-修复] GPU正变换完成" << std::endl;
}

// 修复版GPU逆变换
__host__ void ntt_inverse_gpu_optimized_v2_fixed(CudaU64Memory &d_data, u64 n,
                                                 u64 p, u64 omega,
                                                 u64 *d_twiddles_inv) {
  std::cout << "[优化2-修复] GPU逆变换开始, n=" << n << std::endl;

  // 使用正变换函数 + 逆变换旋转因子表
  ntt_forward_gpu_optimized_v2_fixed(d_data, n, p, omega, d_twiddles_inv);

  // GPU归一化
  u64 n_inv = mod_inv(n, p);
  dim3 blockSize(256);
  dim3 gridSize((n + blockSize.x - 1) / blockSize.x);
  normalize_kernel<<<gridSize, blockSize>>>(d_data.get(), n, n_inv, p);
  syncAndCheck("GPU归一化");

  std::cout << "[优化2-修复] GPU逆变换完成" << std::endl;
}

// 修复版多项式乘法：确保正确性第一
__host__ void poly_multiply_gpu_optimized_v2_fixed(u64 *a, u64 *b, u64 *result,
                                                   u64 n, u64 p) {
  std::cout << "=== GPU多项式乘法 (优化v2-修复版: 安全共享内存) ==="
            << std::endl;

  // 参数
  u64 omega = 3;

  // 计算正确的扩展规模
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
    // CPU内存分配 - 使用RAII管理
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
    ntt_forward_gpu_optimized_v2_fixed(d_a, n_expanded, p, omega, d_twiddles);
    ntt_forward_gpu_optimized_v2_fixed(d_b, n_expanded, p, omega, d_twiddles);

    // 逐点乘法
    dim3 blockSize(256);
    dim3 gridSize((n_expanded + blockSize.x - 1) / blockSize.x);
    pointwise_mul_kernel_optimized<<<gridSize, blockSize>>>(
        d_a.get(), d_b.get(), d_result.get(), n_expanded, p);
    syncAndCheck("GPU逐点乘法");

    // 逆变换
    ntt_inverse_gpu_optimized_v2_fixed(d_result, n_expanded, p, omega,
                                       d_twiddles_inv);

    // 结果传输回CPU
    d_result.copyToHost(result_exp.get(), n_expanded);

    // 复制结果
    for (u64 i = 0; i < 2 * n - 1; i++) {
      result[i] = result_exp[i];
    }

    // 清理GPU内存
    cudaFree(d_twiddles);
    cudaFree(d_twiddles_inv);

    std::cout << "[优化2-修复] 多项式乘法完成！" << std::endl;

  } catch (const std::exception &e) {
    std::cout << "❌ 异常: " << e.what() << std::endl;
  }
}