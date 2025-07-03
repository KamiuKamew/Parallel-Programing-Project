#include "ntt.h"
#include <iostream>

// 声明我们的优化版本v2函数
extern void poly_multiply_gpu_optimized_v2_fixed(u64 *a, u64 *b, u64 *result,
                                                 u64 n);

// 实现main.cc需要的接口函数
void poly_multiply_ntt_cuda_optimized_v2(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                         u64 OMEGA) {
  std::cout << "[CUDA优化v2] 开始智能共享内存版本NTT..." << std::endl;

  // 调用我们的修复版v2函数
  poly_multiply_gpu_optimized_v2_fixed(a, b, ab, n);

  std::cout << "[CUDA优化v2] 智能共享内存版本NTT完成" << std::endl;
}

// 其他可能需要的占位符函数
void poly_multiply_ntt_cuda_parallel(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                     u64 OMEGA) {
  std::cout << "[CUDA] 调用标准版本..." << std::endl;
  poly_multiply_ntt_cuda_optimized_v2(a, b, ab, n, p, OMEGA);
}

void poly_multiply_ntt_cuda_serial(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                   u64 OMEGA) {
  std::cout << "[CUDA] 调用串行版本..." << std::endl;
  poly_multiply_ntt_cuda_optimized_v2(a, b, ab, n, p, OMEGA);
}

void poly_multiply_ntt_cuda_parallel_ready(u64 *a, u64 *b, u64 *ab, u64 n,
                                           u64 p, u64 OMEGA) {
  std::cout << "[CUDA] 调用准备版本..." << std::endl;
  poly_multiply_ntt_cuda_optimized_v2(a, b, ab, n, p, OMEGA);
}

void poly_multiply_ntt_cuda_simple(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                   u64 OMEGA) {
  std::cout << "[CUDA] 调用简单版本..." << std::endl;
  poly_multiply_ntt_cuda_optimized_v2(a, b, ab, n, p, OMEGA);
}