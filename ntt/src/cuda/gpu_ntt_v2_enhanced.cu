#include "gpu_ntt_common.h"
#include <cuda_runtime.h>
#include <iostream>

// 包含v1版本的实现作为基础
#include "gpu_ntt_v1.cu"

// 性能增强版v2：优化线程配置和内存访问
#define WARP_SIZE 32
#define OPTIMAL_BLOCK_SIZE 512   // 更大的block size提高占用率
#define SHARED_MEM_ELEMENTS 1024 // 优化共享内存大小

// 增强版蝶形运算kernel：优化内存合并访问
__global__ void butterfly_kernel_enhanced(u64 *data, u64 n, u64 mid,
                                          u64 *twiddles, u64 twiddle_offset,
                                          u64 p) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int total_operations = n / 2;

  if (idx >= total_operations)
    return;

  // 计算(j,k)坐标 - 与v1版本一致确保正确性
  int block_id = idx / mid;
  int k = idx % mid;
  int j = block_id * (mid << 1);

  int i1 = j + k;
  int i2 = j + k + mid;

  // 预取旋转因子到寄存器中
  u64 w = twiddles[twiddle_offset + k];

  // 预取数据到寄存器中减少全局内存访问
  u64 u_val = data[i1];
  u64 v_val = data[i2];

  // 蝶形运算
  u64 v_mul = mod_mul(v_val, w, p);
  data[i1] = mod_add(u_val, v_mul, p);
  data[i2] = mod_sub(u_val, v_mul, p);
}

// 增强版共享内存kernel：针对小规模mid优化
__global__ void butterfly_kernel_shared_enhanced(u64 *data, u64 n, u64 mid,
                                                 u64 *twiddles,
                                                 u64 twiddle_offset, u64 p) {
  extern __shared__ u64 shared_mem[];

  int tid = threadIdx.x;
  int bid = blockIdx.x;
  int block_size = blockDim.x;

  // 计算每个block处理的蝶形运算数量
  int operations_per_block = SHARED_MEM_ELEMENTS / 2;
  int start_op = bid * operations_per_block;
  int end_op = min(start_op + operations_per_block, (int)(n / 2));

  if (start_op >= n / 2)
    return;

  // 协作加载数据到共享内存
  for (int op = start_op + tid; op < end_op; op += block_size) {
    int block_id = op / mid;
    int k = op % mid;
    int j = block_id * (mid << 1);

    int i1 = j + k;
    int i2 = j + k + mid;

    // 加载到共享内存（使用局部索引）
    int local_i1 = (i1 - start_op * 2) % SHARED_MEM_ELEMENTS;
    int local_i2 = (i2 - start_op * 2) % SHARED_MEM_ELEMENTS;

    if (local_i1 < SHARED_MEM_ELEMENTS && local_i2 < SHARED_MEM_ELEMENTS) {
      shared_mem[local_i1] = data[i1];
      shared_mem[local_i2] = data[i2];
    }
  }

  __syncthreads();

  // 在共享内存中执行蝶形运算
  for (int op = start_op + tid; op < end_op; op += block_size) {
    int block_id = op / mid;
    int k = op % mid;
    int j = block_id * (mid << 1);

    int i1 = j + k;
    int i2 = j + k + mid;

    int local_i1 = (i1 - start_op * 2) % SHARED_MEM_ELEMENTS;
    int local_i2 = (i2 - start_op * 2) % SHARED_MEM_ELEMENTS;

    if (local_i1 < SHARED_MEM_ELEMENTS && local_i2 < SHARED_MEM_ELEMENTS) {
      u64 w = twiddles[twiddle_offset + k];

      u64 u_val = shared_mem[local_i1];
      u64 v_val = shared_mem[local_i2];
      u64 v_mul = mod_mul(v_val, w, p);

      shared_mem[local_i1] = mod_add(u_val, v_mul, p);
      shared_mem[local_i2] = mod_sub(u_val, v_mul, p);
    }
  }

  __syncthreads();

  // 写回到全局内存
  for (int op = start_op + tid; op < end_op; op += block_size) {
    int block_id = op / mid;
    int k = op % mid;
    int j = block_id * (mid << 1);

    int i1 = j + k;
    int i2 = j + k + mid;

    int local_i1 = (i1 - start_op * 2) % SHARED_MEM_ELEMENTS;
    int local_i2 = (i2 - start_op * 2) % SHARED_MEM_ELEMENTS;

    if (local_i1 < SHARED_MEM_ELEMENTS && local_i2 < SHARED_MEM_ELEMENTS) {
      data[i1] = shared_mem[local_i1];
      data[i2] = shared_mem[local_i2];
    }
  }
}

// 增强版GPU正变换：智能选择kernel策略
__host__ void ntt_forward_gpu_enhanced(CudaU64Memory &d_data, u64 n, u64 p,
                                       u64 omega, u64 *d_twiddles) {
  std::cout << "[增强版] GPU正变换开始, n=" << n << std::endl;

  // 位逆序排列 - 使用优化的线程配置
  int bits = 0;
  u64 temp = n - 1;
  while (temp) {
    bits++;
    temp >>= 1;
  }

  dim3 blockSize(OPTIMAL_BLOCK_SIZE);
  dim3 gridSize((n + blockSize.x - 1) / blockSize.x);
  bit_reverse_kernel_optimized<<<gridSize, blockSize>>>(d_data.get(), n, bits);
  syncAndCheck("GPU位逆序排列");

  // 蝶形运算 - 智能选择kernel
  u64 twiddle_offset = 0;
  for (u64 mid = 1; mid < n; mid <<= 1) {
    int total_operations = n / 2;

    // 根据mid大小智能选择kernel
    if (mid <= 64 && n <= 65536) {
      // 小规模mid：使用共享内存优化
      int operations_per_block = SHARED_MEM_ELEMENTS / 2;
      dim3 sharedGrid((total_operations + operations_per_block - 1) /
                      operations_per_block);
      dim3 sharedBlock((OPTIMAL_BLOCK_SIZE < operations_per_block)
                           ? OPTIMAL_BLOCK_SIZE
                           : operations_per_block);

      butterfly_kernel_shared_enhanced<<<sharedGrid, sharedBlock,
                                         SHARED_MEM_ELEMENTS * sizeof(u64)>>>(
          d_data.get(), n, mid, d_twiddles, twiddle_offset, p);
      syncAndCheck("增强版共享内存蝶形运算");
    } else {
      // 大规模mid：使用全局内存优化版本
      dim3 enhancedGrid((total_operations + OPTIMAL_BLOCK_SIZE - 1) /
                        OPTIMAL_BLOCK_SIZE);
      dim3 enhancedBlock(OPTIMAL_BLOCK_SIZE);

      butterfly_kernel_enhanced<<<enhancedGrid, enhancedBlock>>>(
          d_data.get(), n, mid, d_twiddles, twiddle_offset, p);
      syncAndCheck("增强版全局内存蝶形运算");
    }

    twiddle_offset += mid;
  }

  std::cout << "[增强版] GPU正变换完成" << std::endl;
}

// 增强版GPU逆变换
__host__ void ntt_inverse_gpu_enhanced(CudaU64Memory &d_data, u64 n, u64 p,
                                       u64 omega, u64 *d_twiddles_inv) {
  std::cout << "[增强版] GPU逆变换开始, n=" << n << std::endl;

  // 使用正变换函数 + 逆变换旋转因子表
  ntt_forward_gpu_enhanced(d_data, n, p, omega, d_twiddles_inv);

  // GPU归一化 - 使用优化的线程配置
  u64 n_inv = mod_inv(n, p);
  dim3 blockSize(OPTIMAL_BLOCK_SIZE);
  dim3 gridSize((n + blockSize.x - 1) / blockSize.x);
  normalize_kernel<<<gridSize, blockSize>>>(d_data.get(), n, n_inv, p);
  syncAndCheck("GPU归一化");

  std::cout << "[增强版] GPU逆变换完成" << std::endl;
}

// 增强版多项式乘法：性能优化版本
__host__ void poly_multiply_gpu_optimized_v2_enhanced(u64 *a, u64 *b,
                                                      u64 *result, u64 n,
                                                      u64 p) {
  std::cout << "=== GPU多项式乘法 (增强版v2: 性能优化) ===" << std::endl;

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
    ntt_forward_gpu_enhanced(d_a, n_expanded, p, omega, d_twiddles);
    ntt_forward_gpu_enhanced(d_b, n_expanded, p, omega, d_twiddles);

    // 逐点乘法 - 使用优化的线程配置
    dim3 blockSize(OPTIMAL_BLOCK_SIZE);
    dim3 gridSize((n_expanded + blockSize.x - 1) / blockSize.x);
    pointwise_mul_kernel_optimized<<<gridSize, blockSize>>>(
        d_a.get(), d_b.get(), d_result.get(), n_expanded, p);
    syncAndCheck("GPU逐点乘法");

    // 逆变换
    ntt_inverse_gpu_enhanced(d_result, n_expanded, p, omega, d_twiddles_inv);

    // 结果传输回CPU
    d_result.copyToHost(result_exp.get(), n_expanded);

    // 复制结果
    for (u64 i = 0; i < 2 * n - 1; i++) {
      result[i] = result_exp[i];
    }

    // 清理GPU内存
    cudaFree(d_twiddles);
    cudaFree(d_twiddles_inv);

    std::cout << "[增强版] 多项式乘法完成！" << std::endl;

  } catch (const std::exception &e) {
    std::cout << "❌ 异常: " << e.what() << std::endl;
  }
}