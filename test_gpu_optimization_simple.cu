#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

// 直接包含必要的头文件
#include "ntt/src/include/general/op.h"
#include "ntt/src/include/general/type.h"
#include "ntt/src/include/general/utils.h"
#include "ntt/src/include/transform.h"

using namespace std::chrono;

// CPU版本NTT实现（简化版本）
template <typename T>
void simple_ntt_cpu(T *a, T *b, T *ab, T n, T p, T omega = 3) {
  using T_mont = T;
  MontMod<T> montMod(p);

  T n_expanded = expand_n(2 * n - 1);
  T *a_expanded = expand_a((T *)a, n, n_expanded);
  T *b_expanded = expand_a((T *)b, n, n_expanded);

  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  T_mont *a_mont = new T_mont[n_expanded]{};
  T_mont *b_mont = new T_mont[n_expanded]{};
  T_mont *ab_mont = new T_mont[n_expanded]{};
  for (T i = 0; i < n_expanded; ++i)
    a_mont[i] = montMod.from_T(a_expanded[i]);
  for (T i = 0; i < n_expanded; ++i)
    b_mont[i] = montMod.from_T(b_expanded[i]);
  T_mont omega_mont = montMod.from_T(omega);

  ntt_forward_mont(a_mont, n_expanded, p, omega_mont);
  ntt_forward_mont(b_mont, n_expanded, p, omega_mont);

  for (T i = 0; i < n_expanded; ++i)
    ab_mont[i] = montMod.mul(a_mont[i], b_mont[i]);

  ntt_inverse_mont(ab_mont, n_expanded, p, montMod.inv(omega_mont));

  for (T i = 0; i < n_expanded; ++i)
    ab[i] = montMod.to_T(ab_mont[i]);

  bit_reverse_permute((T *)ab, n_expanded);

  delete[] a_expanded;
  delete[] b_expanded;
  delete[] a_mont;
  delete[] b_mont;
  delete[] ab_mont;
}

// GPU版本函数声明
extern void poly_multiply_ntt_gpu_mont_optimized(uint32_t *a, uint32_t *b,
                                                 uint32_t *ab, uint32_t n,
                                                 uint32_t p, uint32_t omega);

void test_optimization_effect() {
  const uint32_t p = 998244353; // NTT模数
  const uint32_t omega = 3;     // 原根

  std::vector<uint32_t> test_sizes = {1024, 4096, 16384, 65536, 131072};

  std::cout << "GPU NTT优化效果测试" << std::endl;
  std::cout << "==================" << std::endl;
  std::cout << std::setw(10) << "Size" << std::setw(15) << "CPU Time(μs)"
            << std::setw(15) << "GPU Time(μs)" << std::setw(15) << "Speedup"
            << std::setw(15) << "Status" << std::endl;
  std::cout << std::string(75, '-') << std::endl;

  for (uint32_t n : test_sizes) {
    // 分配测试数据
    uint32_t *a = new uint32_t[n];
    uint32_t *b = new uint32_t[n];
    uint32_t *cpu_result = new uint32_t[2 * n - 1];
    uint32_t *gpu_result = new uint32_t[2 * n - 1];

    // 初始化为0
    for (uint32_t i = 0; i < 2 * n - 1; i++) {
      cpu_result[i] = 0;
      gpu_result[i] = 0;
    }

    // 生成测试数据
    for (uint32_t i = 0; i < n; i++) {
      a[i] = (i + 1) % 1000; // 使用小一点的数值避免溢出
      b[i] = (i * 2 + 1) % 1000;
    }

    std::cout << "Testing n=" << n << "..." << std::flush;

    // CPU测试
    auto cpu_start = high_resolution_clock::now();
    simple_ntt_cpu(a, b, cpu_result, n, p, omega);
    auto cpu_end = high_resolution_clock::now();
    double cpu_time = duration_cast<microseconds>(cpu_end - cpu_start).count();

    // GPU优化版本测试
    auto gpu_start = high_resolution_clock::now();
    poly_multiply_ntt_gpu_mont_optimized(a, b, gpu_result, n, p, omega);
    auto gpu_end = high_resolution_clock::now();
    double gpu_time = duration_cast<microseconds>(gpu_end - gpu_start).count();

    // 验证正确性
    bool correct = true;
    for (uint32_t i = 0; i < 2 * n - 1; i++) {
      if (cpu_result[i] != gpu_result[i]) {
        correct = false;
        std::cout << "\rMismatch at index " << i << ": CPU=" << cpu_result[i]
                  << ", GPU=" << gpu_result[i] << std::endl;
        break;
      }
    }

    double speedup = cpu_time / gpu_time;

    std::cout << "\r" << std::setw(10) << n << std::setw(15) << std::fixed
              << std::setprecision(1) << cpu_time << std::setw(15) << std::fixed
              << std::setprecision(1) << gpu_time << std::setw(15) << std::fixed
              << std::setprecision(3) << speedup << std::setw(15)
              << (correct ? "PASS" : "FAIL");
    std::cout << std::endl;

    delete[] a;
    delete[] b;
    delete[] cpu_result;
    delete[] gpu_result;
  }
}

int main() {
  // 检查CUDA是否可用
  int deviceCount;
  cudaError_t error = cudaGetDeviceCount(&deviceCount);

  if (error != cudaSuccess || deviceCount == 0) {
    std::cout << "No CUDA device found!" << std::endl;
    return 1;
  }

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  std::cout << "Using GPU: " << prop.name << std::endl;
  std::cout << "Compute Capability: " << prop.major << "." << prop.minor
            << std::endl;
  std::cout << std::endl;

  test_optimization_effect();
  return 0;
}