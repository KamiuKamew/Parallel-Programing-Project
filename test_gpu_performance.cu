#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

// 包含必要的头文件
#include "ntt/src/include/general/op.h"
#include "ntt/src/include/general/type.h"
#include "ntt/src/include/general/utils.h"
#include "ntt/src/include/transform.h"

using namespace std::chrono;

// 先包含GPU实现代码，这样就有了所有函数声明
#include "ntt/src/include/CUDA/ntt_gpu_optimized.cu"

// CPU版本NTT实现（完整版本）
template <typename T>
void cpu_ntt_baseline(T *a, T *b, T *ab, T n, T p, T omega = 3) {
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

struct TestResult {
  uint32_t n;
  double cpu_time_us;
  double gpu_time_us;
  double speedup;
  bool correct;
};

void run_performance_test() {
  const uint32_t p = 998244353; // NTT模数
  const uint32_t omega = 3;     // 原根

  // 测试不同规模
  std::vector<uint32_t> test_sizes = {1024, 4096, 16384, 65536, 131072, 262144};
  std::vector<TestResult> results;

  std::cout << "GPU NTT 优化效果验证测试" << std::endl;
  std::cout << "========================" << std::endl;
  std::cout << std::setw(10) << "Size" << std::setw(15) << "CPU Time(μs)"
            << std::setw(15) << "GPU Time(μs)" << std::setw(12) << "Speedup"
            << std::setw(10) << "Status" << std::endl;
  std::cout << std::string(62, '-') << std::endl;

  for (uint32_t n : test_sizes) {
    TestResult result;
    result.n = n;

    // 分配内存
    uint32_t *a = new uint32_t[n];
    uint32_t *b = new uint32_t[n];
    uint32_t *cpu_result = new uint32_t[2 * n - 1];
    uint32_t *gpu_result = new uint32_t[2 * n - 1];

    // 初始化测试数据
    for (uint32_t i = 0; i < n; i++) {
      a[i] = (i + 1) % 1000;
      b[i] = (i * 2 + 3) % 1000;
    }

    // 初始化结果数组
    for (uint32_t i = 0; i < 2 * n - 1; i++) {
      cpu_result[i] = 0;
      gpu_result[i] = 0;
    }

    // 预热GPU
    poly_multiply_ntt_gpu_mont_optimized(a, b, gpu_result, n, p, omega);

    // CPU性能测试
    auto cpu_start = high_resolution_clock::now();
    cpu_ntt_baseline(a, b, cpu_result, n, p, omega);
    auto cpu_end = high_resolution_clock::now();
    result.cpu_time_us =
        duration_cast<microseconds>(cpu_end - cpu_start).count();

    // GPU性能测试
    auto gpu_start = high_resolution_clock::now();
    poly_multiply_ntt_gpu_mont_optimized(a, b, gpu_result, n, p, omega);
    auto gpu_end = high_resolution_clock::now();
    result.gpu_time_us =
        duration_cast<microseconds>(gpu_end - gpu_start).count();

    // 验证正确性
    result.correct = true;
    for (uint32_t i = 0; i < 2 * n - 1; i++) {
      if (cpu_result[i] != gpu_result[i]) {
        result.correct = false;
        break;
      }
    }

    result.speedup = result.cpu_time_us / result.gpu_time_us;
    results.push_back(result);

    // 输出结果
    std::cout << std::setw(10) << result.n << std::setw(15) << std::fixed
              << std::setprecision(1) << result.cpu_time_us << std::setw(15)
              << std::fixed << std::setprecision(1) << result.gpu_time_us
              << std::setw(12) << std::fixed << std::setprecision(3)
              << result.speedup << std::setw(10)
              << (result.correct ? "PASS" : "FAIL");
    std::cout << std::endl;

    delete[] a;
    delete[] b;
    delete[] cpu_result;
    delete[] gpu_result;
  }

  // 输出统计信息
  std::cout << std::endl;
  std::cout << "=== 优化效果统计 ===" << std::endl;

  // 计算各规模的平均加速比
  double total_speedup = 0;
  int valid_tests = 0;
  double max_speedup = 0;
  uint32_t max_speedup_size = 0;

  for (const auto &result : results) {
    if (result.correct) {
      total_speedup += result.speedup;
      valid_tests++;
      if (result.speedup > max_speedup) {
        max_speedup = result.speedup;
        max_speedup_size = result.n;
      }
    }
  }

  if (valid_tests > 0) {
    double avg_speedup = total_speedup / valid_tests;
    std::cout << "平均加速比: " << std::fixed << std::setprecision(3)
              << avg_speedup << "x" << std::endl;
    std::cout << "最大加速比: " << std::fixed << std::setprecision(3)
              << max_speedup << "x (n=" << max_speedup_size << ")" << std::endl;

    // 分析大规模问题的性能
    double large_scale_speedup = 0;
    int large_count = 0;
    for (const auto &result : results) {
      if (result.n >= 65536 && result.correct) {
        large_scale_speedup += result.speedup;
        large_count++;
      }
    }

    if (large_count > 0) {
      large_scale_speedup /= large_count;
      std::cout << "大规模问题平均加速比 (n≥65536): " << std::fixed
                << std::setprecision(3) << large_scale_speedup << "x"
                << std::endl;
    }
  }

  // 生成CSV数据
  std::cout << std::endl;
  std::cout << "CSV格式数据:" << std::endl;
  std::cout << "Size,CPU_Time_us,GPU_Time_us,Speedup,Correct" << std::endl;
  for (const auto &result : results) {
    std::cout << result.n << "," << result.cpu_time_us << ","
              << result.gpu_time_us << "," << result.speedup << ","
              << (result.correct ? "true" : "false") << std::endl;
  }
}

int main() {
  // 检查CUDA设备
  int deviceCount;
  cudaError_t error = cudaGetDeviceCount(&deviceCount);

  if (error != cudaSuccess || deviceCount == 0) {
    std::cerr << "错误: 未找到CUDA设备!" << std::endl;
    return 1;
  }

  // 输出GPU信息
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  std::cout << "GPU设备: " << prop.name << std::endl;
  std::cout << "计算能力: " << prop.major << "." << prop.minor << std::endl;
  std::cout << "SM数量: " << prop.multiProcessorCount << std::endl;
  std::cout << "全局内存: " << prop.totalGlobalMem / (1024 * 1024) << " MB"
            << std::endl;
  std::cout << std::endl;

  try {
    run_performance_test();
  } catch (const std::exception &e) {
    std::cerr << "测试过程中发生错误: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}

// 包含GPU实现代码
#include "ntt/src/include/CUDA/ntt_gpu_optimized.cu"