#include "ntt/src/include/ntt.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std::chrono;

void test_optimization_effect() {
  const uint32_t p = 998244353; // NTT模数
  const uint32_t omega = 3;     // 原根

  std::vector<uint32_t> test_sizes = {1024, 4096, 16384, 65536, 131072, 262144};

  std::cout << "GPU NTT优化效果测试" << std::endl;
  std::cout << "==================" << std::endl;
  std::cout << std::setw(10) << "Size" << std::setw(15) << "CPU Time(μs)"
            << std::setw(15) << "GPU Time(μs)" << std::setw(15) << "Speedup"
            << std::endl;
  std::cout << std::string(60, '-') << std::endl;

  for (uint32_t n : test_sizes) {
    // 分配测试数据
    uint32_t *a = new uint32_t[n];
    uint32_t *b = new uint32_t[n];
    uint32_t *cpu_result = new uint32_t[2 * n - 1];
    uint32_t *gpu_result = new uint32_t[2 * n - 1];

    // 生成测试数据
    for (uint32_t i = 0; i < n; i++) {
      a[i] = (i + 1) % p;
      b[i] = (i * 2 + 1) % p;
    }

    // CPU测试
    auto cpu_start = high_resolution_clock::now();
    poly_multiply_ntt(a, b, cpu_result, n, p, omega);
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
        break;
      }
    }

    double speedup = cpu_time / gpu_time;

    std::cout << std::setw(10) << n << std::setw(15) << std::fixed
              << std::setprecision(1) << cpu_time << std::setw(15) << std::fixed
              << std::setprecision(1) << gpu_time << std::setw(15) << std::fixed
              << std::setprecision(3) << speedup;

    if (!correct) {
      std::cout << " [INCORRECT]";
    }
    std::cout << std::endl;

    delete[] a;
    delete[] b;
    delete[] cpu_result;
    delete[] gpu_result;
  }
}

int main() {
  test_optimization_effect();
  return 0;
}