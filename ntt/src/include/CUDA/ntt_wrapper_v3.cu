#include "ntt.h"
#include <iostream>

extern void poly_multiply_gpu_optimized_v3_test3(u64 *a, u64 *b, u64 *result,
                                                 u64 n);

void poly_multiply_ntt_cuda_optimized_v3(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                         u64 OMEGA = 3) {
  std::cout << "[CUDA优化v3] 调用多层共享内存实现..." << std::endl;
  poly_multiply_gpu_optimized_v3_test3(a, b, ab, n);
  std::cout << "[CUDA优化v3] 完成" << std::endl;
}