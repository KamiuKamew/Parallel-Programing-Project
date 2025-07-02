#include "ntt.h"
#include <cassert>
#include <iostream>

// CUDA错误检查辅助函数（直接使用u64）
void check_cuda_error_simple(const char *msg) {
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::cerr << "CUDA Error (" << msg << "): " << cudaGetErrorString(err)
              << std::endl;
    exit(EXIT_FAILURE);
  }
}

// 简化的GPU串行NTT - 第一步：验证CUDA环境（直接使用u64）
__global__ void simple_ntt_kernel(u64 *a, u64 *b, u64 *result, u64 n, u64 p) {
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    // 最简化的实现：直接进行朴素多项式乘法验证GPU运算能力
    for (u64 i = 0; i < n + n - 1; ++i) {
      result[i] = 0;
      for (u64 j = 0; j < n; ++j) {
        for (u64 k = 0; k < n; ++k) {
          if (j + k == i) {
            result[i] = (result[i] + (a[j] * b[k]) % p) % p;
          }
        }
      }
    }
  }
}

// 第一步：最简化的CUDA版本 - 验证GPU环境（直接使用u64）
void poly_multiply_ntt_cuda_simple(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                   u64 OMEGA) {
  std::cout << "[CUDA Simple] 开始最简化版本..." << std::endl;

  // 分配GPU内存
  u64 *d_a, *d_b, *d_result;
  size_t result_size = 2 * n - 1;

  cudaMalloc((void **)&d_a, n * sizeof(u64));
  cudaMalloc((void **)&d_b, n * sizeof(u64));
  cudaMalloc((void **)&d_result, result_size * sizeof(u64));
  check_cuda_error_simple("GPU memory allocation");

  // 复制数据到GPU
  cudaMemcpy(d_a, a, n * sizeof(u64), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, b, n * sizeof(u64), cudaMemcpyHostToDevice);
  check_cuda_error_simple("copy to GPU");

  // 执行简化的GPU计算
  simple_ntt_kernel<<<1, 1>>>(d_a, d_b, d_result, n, p);
  cudaDeviceSynchronize();
  check_cuda_error_simple("kernel execution");

  // 复制结果回CPU
  cudaMemcpy(ab, d_result, result_size * sizeof(u64), cudaMemcpyDeviceToHost);
  check_cuda_error_simple("copy to CPU");

  // 清理GPU内存
  cudaFree(d_a);
  cudaFree(d_b);
  cudaFree(d_result);
  check_cuda_error_simple("GPU memory cleanup");

  std::cout << "[CUDA Simple] 最简化版本完成" << std::endl;
}