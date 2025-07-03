#include "gpu_ntt_common.h"
#include <cuda_runtime.h>
#include <iostream>

// 包含v1版本的实现
#include "gpu_ntt_v1.cu"

// 共享内存优化版本：将多个蝶形层合并到共享内存中处理
#define SHARED_MEM_SIZE 1024 // 共享内存大小限制
#define MAX_THREADS_PER_BLOCK 512

// 声明需要的函数 (从v1版本引用，避免重复定义)
extern void precompute_twiddle_factors(u64 **d_twiddles, u64 n, u64 p,
                                       u64 omega);
extern __global__ void bit_reverse_kernel_optimized(u64 *data, u64 n, int bits);
extern __global__ void butterfly_kernel_optimized(u64 *data, u64 n, u64 mid,
                                                  u64 *twiddles,
                                                  u64 twiddle_offset, u64 p);
extern __global__ void normalize_kernel(u64 *data, u64 n, u64 n_inv, u64 p);
extern __global__ void
pointwise_mul_kernel_optimized(u64 *a, u64 *b, u64 *result, u64 n, u64 p);

// 新的共享内存优化kernel
__global__ void butterfly_kernel_shared_memory(u64 *data, u64 n, u64 start_mid,
                                               u64 end_mid, u64 *twiddles,
                                               u64 twiddle_offset, u64 p) {
  // 计算共享内存大小
  extern __shared__ u64 shared_data[];

  int tid = threadIdx.x;
  int bid = blockIdx.x;
  int block_size = blockDim.x;

  // 每个block处理的数据范围
  int elements_per_block = SHARED_MEM_SIZE;
  int start_idx = bid * elements_per_block;

  if (start_idx >= n)
    return;

  int end_idx = min(start_idx + elements_per_block, (int)n);
  int actual_elements = end_idx - start_idx;

  // 将数据从全局内存复制到共享内存
  for (int i = tid; i < actual_elements; i += block_size) {
    if (start_idx + i < n) {
      shared_data[i] = data[start_idx + i];
    }
  }
  __syncthreads();

  // 在共享内存中执行多层蝶形运算
  u64 current_twiddle_offset = twiddle_offset;

  for (u64 mid = start_mid; mid <= end_mid; mid <<= 1) {
    // 确保mid层级在当前block的处理范围内
    if (mid > actual_elements / 2)
      break;

    for (int stride = tid; stride < actual_elements / 2; stride += block_size) {
      int block_id = stride / mid;
      int k = stride % mid;
      int j = block_id * (mid << 1);

      int i1 = j + k;
      int i2 = j + k + mid;

      // 检查索引边界
      if (i2 < actual_elements) {
        u64 w = twiddles[current_twiddle_offset + k];

        u64 u_val = shared_data[i1];
        u64 v_val = mod_mul(shared_data[i2], w, p);
        shared_data[i1] = mod_add(u_val, v_val, p);
        shared_data[i2] = mod_sub(u_val, v_val, p);
      }
    }

    __syncthreads();
    current_twiddle_offset += mid;
  }

  // 将结果从共享内存复制回全局内存
  for (int i = tid; i < actual_elements; i += block_size) {
    if (start_idx + i < n) {
      data[start_idx + i] = shared_data[i];
    }
  }
}

// 优化版本v2：使用共享内存的GPU正变换
__host__ void ntt_forward_gpu_optimized_v2(CudaU64Memory &d_data, u64 n, u64 p,
                                           u64 omega, u64 *d_twiddles) {
  std::cout << "[优化2] GPU正变换开始 (共享内存), n=" << n << std::endl;

  // 位逆序排列
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

  // 使用共享内存优化的蝶形运算
  u64 twiddle_offset = 0;

  // 分组处理：将小的mid层级合并到共享内存中处理
  for (u64 group_start = 1; group_start < n;) {
    u64 group_end = group_start;
    u64 total_shared_elements = 0;

    // 计算可以在共享内存中处理的层级范围
    while (group_end < n && total_shared_elements < SHARED_MEM_SIZE / 4) {
      total_shared_elements += group_end;
      group_end <<= 1;
    }

    if (group_end > n)
      group_end = n;

    // 对于可以在共享内存中处理的层级，使用共享内存kernel
    if (total_shared_elements < SHARED_MEM_SIZE / 4) {
      // 计算共享内存大小
      size_t shared_mem_bytes = SHARED_MEM_SIZE * sizeof(u64);

      // 计算grid尺寸
      int elements_per_block = SHARED_MEM_SIZE;
      dim3 sharedGrid((n + elements_per_block - 1) / elements_per_block);
      dim3 sharedBlock(256);

      butterfly_kernel_shared_memory<<<sharedGrid, sharedBlock,
                                       shared_mem_bytes>>>(
          d_data.get(), n, group_start, group_end >> 1, d_twiddles,
          twiddle_offset, p);
      syncAndCheck("共享内存蝶形运算");

      // 更新偏移量
      for (u64 mid = group_start; mid < group_end; mid <<= 1) {
        twiddle_offset += mid;
      }
    } else {
      // 对于大的层级，回退到原始kernel
      u64 mid = group_start;
      dim3 blockSize(256);
      dim3 gridSize((n / 2 + blockSize.x - 1) / blockSize.x);
      butterfly_kernel_optimized<<<gridSize, blockSize>>>(
          d_data.get(), n, mid, d_twiddles, twiddle_offset, p);
      syncAndCheck("标准蝶形运算");
      twiddle_offset += mid;
    }

    group_start = group_end;
  }

  std::cout << "[优化2] GPU正变换完成" << std::endl;
}

// 优化版本v2：使用共享内存的GPU逆变换
__host__ void ntt_inverse_gpu_optimized_v2(CudaU64Memory &d_data, u64 n, u64 p,
                                           u64 omega, u64 *d_twiddles_inv) {
  std::cout << "[优化2] GPU逆变换开始 (共享内存), n=" << n << std::endl;

  // 使用正变换函数 + 逆变换旋转因子表
  ntt_forward_gpu_optimized_v2(d_data, n, p, omega, d_twiddles_inv);

  // GPU归一化
  u64 n_inv = mod_inv(n, p);
  dim3 blockSize(256);
  dim3 gridSize((n + blockSize.x - 1) / blockSize.x);
  normalize_kernel<<<gridSize, blockSize>>>(d_data.get(), n, n_inv, p);
  syncAndCheck("GPU归一化");

  std::cout << "[优化2] GPU逆变换完成" << std::endl;
}

// 优化版本v2：完整的多项式乘法
__host__ void poly_multiply_gpu_optimized_v2(u64 *a, u64 *b, u64 *result, u64 n,
                                             u64 p) {
  std::cout << "=== GPU多项式乘法 (优化v2: 共享内存) ===" << std::endl;

  // 参数 - omega固定为3（题目保证所有模数的原根都是3）
  u64 omega = 3;

  // 计算正确的扩展规模：找到满足 n_expanded >= 2*n-1 的最小2的幂次
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

  // 分配GPU内存
  CudaU64Memory d_a(n_expanded);
  CudaU64Memory d_b(n_expanded);
  CudaU64Memory d_result(n_expanded);

  // 复制数据到GPU
  cudaMemcpy(d_a.get(), a, n * sizeof(u64), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b.get(), b, n * sizeof(u64), cudaMemcpyHostToDevice);
  cudaMemset(d_a.get() + n, 0, n * sizeof(u64));
  cudaMemset(d_b.get() + n, 0, n * sizeof(u64));

  // 预计算旋转因子
  u64 *d_twiddles_forward, *d_twiddles_inverse;
  precompute_twiddle_factors(&d_twiddles_forward, n_expanded, p, omega);
  precompute_twiddle_factors(&d_twiddles_inverse, n_expanded, p,
                             mod_inv(omega, p));

  // 正变换
  ntt_forward_gpu_optimized_v2(d_a, n_expanded, p, omega, d_twiddles_forward);
  ntt_forward_gpu_optimized_v2(d_b, n_expanded, p, omega, d_twiddles_forward);

  // 逐点乘法
  dim3 blockSize(256);
  dim3 gridSize((n_expanded + blockSize.x - 1) / blockSize.x);
  pointwise_mul_kernel_optimized<<<gridSize, blockSize>>>(
      d_a.get(), d_b.get(), d_result.get(), n_expanded, p);
  syncAndCheck("GPU逐点乘法");

  // 逆变换
  ntt_inverse_gpu_optimized_v2(d_result, n_expanded, p, omega,
                               d_twiddles_inverse);

  // 复制结果回CPU
  cudaMemcpy(result, d_result.get(), n_expanded * sizeof(u64),
             cudaMemcpyDeviceToHost);

  // 清理
  cudaFree(d_twiddles_forward);
  cudaFree(d_twiddles_inverse);

  std::cout << "=== GPU多项式乘法完成 (优化v2) ===" << std::endl;
}