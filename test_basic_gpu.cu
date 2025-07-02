#include "ntt/src/include/CUDA/ntt_gpu.cu" // 只包含基础GPU实现
#include "ntt/src/include/ntt.h"
#include <chrono>
#include <iomanip>
#include <iostream>

using namespace std::chrono;

template <typename T> void generate_test_data(T *a, T *b, T n) {
  for (T i = 0; i < n; i++) {
    a[i] = (i * 17 + 7) % 1000 + 1;
    b[i] = (i * 23 + 11) % 1000 + 1;
  }
}

void test_basic_gpu_functionality() {
  const uint32_t p = 998244353;
  const uint32_t omega = 3;

  std::vector<uint32_t> test_sizes = {4096, 16384, 65536};

  std::cout << "=== 基础GPU功能测试 ===" << std::endl;
  std::cout << std::setw(10) << "Size" << std::setw(15) << "CPU时间(μs)"
            << std::setw(15) << "GPU时间(μs)" << std::setw(12) << "加速比"
            << std::setw(10) << "正确性" << std::endl;
  std::cout << std::string(60, '-') << std::endl;

  for (uint32_t n : test_sizes) {
    uint32_t *a = new uint32_t[n];
    uint32_t *b = new uint32_t[n];
    uint32_t *cpu_result = new uint32_t[2 * n - 1]();
    uint32_t *gpu_result = new uint32_t[2 * n - 1]();

    generate_test_data(a, b, n);

    // CPU版本测试
    uint32_t *a_cpu = new uint32_t[n];
    uint32_t *b_cpu = new uint32_t[n];
    std::copy(a, a + n, a_cpu);
    std::copy(b, b + n, b_cpu);

    auto cpu_start = high_resolution_clock::now();
    poly_multiply_ntt(a_cpu, b_cpu, cpu_result, n, p, omega);
    auto cpu_end = high_resolution_clock::now();
    double cpu_time = duration_cast<microseconds>(cpu_end - cpu_start).count();

    delete[] a_cpu;
    delete[] b_cpu;

    // GPU版本测试
    uint32_t *a_gpu = new uint32_t[n];
    uint32_t *b_gpu = new uint32_t[n];
    std::copy(a, a + n, a_gpu);
    std::copy(b, b + n, b_gpu);

    auto gpu_start = high_resolution_clock::now();
    poly_multiply_ntt_gpu_mont(a_gpu, b_gpu, gpu_result, n, p, omega);
    auto gpu_end = high_resolution_clock::now();
    double gpu_time = duration_cast<microseconds>(gpu_end - gpu_start).count();

    delete[] a_gpu;
    delete[] b_gpu;

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
              << std::setprecision(1) << gpu_time << std::setw(12) << std::fixed
              << std::setprecision(2) << speedup << std::setw(10)
              << (correct ? "✓" : "✗") << std::endl;

    delete[] a;
    delete[] b;
    delete[] cpu_result;
    delete[] gpu_result;
  }
}

int main() {
  // 检查CUDA设备
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0) {
    std::cerr << "没有找到CUDA设备!" << std::endl;
    return 1;
  }

  // 显示GPU信息
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  std::cout << "GPU设备: " << prop.name << std::endl;
  std::cout << "计算能力: " << prop.major << "." << prop.minor << std::endl;
  std::cout << std::endl;

  try {
    test_basic_gpu_functionality();
    std::cout << std::endl << "基础GPU测试完成!" << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "测试错误: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}