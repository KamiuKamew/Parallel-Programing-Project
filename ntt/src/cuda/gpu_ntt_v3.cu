#include "gpu_ntt_common.h"
#include <cuda_runtime.h>
#include <iostream>
#include <vector>

// 包含v1版本的基础实现
#include "gpu_ntt_v1.cu"

// v3版本的前向变换实现（简化版）
void ntt_forward_gpu_optimized_v3_test3(CudaU64Memory &d_data, u64 n, u64 p,
                                        u64 omega, u64 *d_twiddles) {
  // 使用v1版本的实现作为基础
  ntt_forward_gpu_optimized(d_data, n, p, omega, d_twiddles);
}

// 逆变换实现
__host__ void ntt_inverse_gpu_optimized_v3_test3(CudaU64Memory &d_data, u64 n,
                                                 u64 p, u64 omega,
                                                 u64 *d_twiddles_inv) {
  ntt_forward_gpu_optimized_v3_test3(d_data, n, p, omega, d_twiddles_inv);
  u64 n_inv = mod_inv(n, p);
  dim3 block(256);
  dim3 grid((n + block.x - 1) / block.x);
  normalize_kernel<<<grid, block>>>(d_data.get(), n, n_inv, p);
  syncAndCheck("v3-test3 归一化");
}

// 多项式乘法包装
void poly_multiply_gpu_optimized_v3_test3(u64 *a, u64 *b, u64 *result, u64 n,
                                          u64 p) {
  const u64 omega = 3;

  // 计算正确的扩展规模
  u64 n_exp = 1;
  while (n_exp < 2 * n - 1) {
    n_exp <<= 1;
  }

  // 验证参数有效性
  if ((p - 1) % n_exp != 0) {
    std::cout << "❌ v3错误: (p-1)=" << (p - 1) << "不能被n_exp=" << n_exp
              << "整除" << std::endl;
    return;
  }

  std::vector<u64> host_a(n_exp, 0), host_b(n_exp, 0);
  for (u64 i = 0; i < n; ++i) {
    host_a[i] = a[i];
    host_b[i] = b[i];
  }

  CudaU64Memory d_a(n_exp), d_b(n_exp), d_res(n_exp);
  d_a.copyFromHost(host_a.data(), n_exp);
  d_b.copyFromHost(host_b.data(), n_exp);

  u64 *d_twiddles_fwd = nullptr, *d_twiddles_inv = nullptr;
  precompute_twiddle_factors(&d_twiddles_fwd, n_exp, p, omega);
  precompute_twiddle_factors(&d_twiddles_inv, n_exp, p, mod_inv(omega, p));

  ntt_forward_gpu_optimized_v3_test3(d_a, n_exp, p, omega, d_twiddles_fwd);
  ntt_forward_gpu_optimized_v3_test3(d_b, n_exp, p, omega, d_twiddles_fwd);

  dim3 block(256);
  dim3 grid((n_exp + block.x - 1) / block.x);
  pointwise_mul_kernel_optimized<<<grid, block>>>(d_a.get(), d_b.get(),
                                                  d_res.get(), n_exp, p);
  syncAndCheck("v3-test3 点乘");

  ntt_inverse_gpu_optimized_v3_test3(d_res, n_exp, p, omega, d_twiddles_inv);

  std::vector<u64> host_res(n_exp);
  d_res.copyToHost(host_res.data(), n_exp);
  for (u64 i = 0; i < n_exp; ++i)
    result[i] = host_res[i];

  cudaFree(d_twiddles_fwd);
  cudaFree(d_twiddles_inv);
}
