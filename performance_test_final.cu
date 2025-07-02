#include "ntt/src/include/ntt.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std::chrono;

// 先包含基础GPU实现
#include "ntt/src/include/CUDA/ntt_gpu.cu"

// 再包含优化GPU实现（去重复定义版本）
#include "ntt/src/include/CUDA/ntt_gpu_optimized.cu"

struct BenchmarkResult {
  uint32_t n;
  double cpu_time_us;
  double gpu_basic_time_us;
  double gpu_opt_time_us;
  double speedup_basic;
  double speedup_opt;
  double improvement_ratio;
  bool correctness_pass;
};

template <typename T> void generate_test_data(T *a, T *b, T n) {
  for (T i = 0; i < n; i++) {
    a[i] = (i * 17 + 7) % 1000 + 1;
    b[i] = (i * 23 + 11) % 1000 + 1;
  }
}

template <typename T>
bool verify_results(T *reference, T *test1, T *test2, T result_size,
                    const std::string &test_name) {
  bool correct1 = true, correct2 = true;

  for (T i = 0; i < result_size; i++) {
    if (reference[i] != test1[i]) {
      correct1 = false;
    }
    if (reference[i] != test2[i]) {
      correct2 = false;
    }
  }

  std::cout << test_name << " - 基础GPU: " << (correct1 ? "✓" : "✗")
            << ", 优化GPU: " << (correct2 ? "✓" : "✗") << std::endl;

  return correct1 && correct2;
}

void run_comprehensive_benchmark() {
  const uint32_t p = 998244353;
  const uint32_t omega = 3;

  std::vector<uint32_t> test_sizes = {1024, 4096, 16384, 65536, 131072, 262144};
  std::vector<BenchmarkResult> results;

  std::cout << "=== GPU NTT 全面性能对比测试 ===" << std::endl;
  std::cout << "对比：CPU基线版本 vs GPU基础版本 vs GPU优化版本" << std::endl;
  std::cout << std::string(85, '=') << std::endl;
  std::cout << std::setw(8) << "Size" << std::setw(12) << "CPU(μs)"
            << std::setw(12) << "GPU基础(μs)" << std::setw(12) << "GPU优化(μs)"
            << std::setw(10) << "基础倍数" << std::setw(10) << "优化倍数"
            << std::setw(10) << "改进比" << std::setw(8) << "正确性"
            << std::endl;
  std::cout << std::string(85, '-') << std::endl;

  for (uint32_t n : test_sizes) {
    BenchmarkResult result;
    result.n = n;

    // 分配内存
    uint32_t *a = new uint32_t[n];
    uint32_t *b = new uint32_t[n];
    uint32_t *cpu_result = new uint32_t[2 * n - 1]();
    uint32_t *gpu_basic_result = new uint32_t[2 * n - 1]();
    uint32_t *gpu_opt_result = new uint32_t[2 * n - 1]();

    // 生成测试数据
    generate_test_data(a, b, n);

    // === CPU版本测试 ===
    const int cpu_runs = (n <= 16384) ? 3 : 1;
    double cpu_total_time = 0.0;

    for (int run = 0; run < cpu_runs; run++) {
      uint32_t *a_copy = new uint32_t[n];
      uint32_t *b_copy = new uint32_t[n];
      std::copy(a, a + n, a_copy);
      std::copy(b, b + n, b_copy);

      auto cpu_start = high_resolution_clock::now();
      poly_multiply_ntt(a_copy, b_copy, cpu_result, n, p, omega);
      auto cpu_end = high_resolution_clock::now();

      cpu_total_time +=
          duration_cast<microseconds>(cpu_end - cpu_start).count();

      delete[] a_copy;
      delete[] b_copy;
    }
    result.cpu_time_us = cpu_total_time / cpu_runs;

    // === GPU基础版本测试 ===
    // 预热
    poly_multiply_ntt_gpu_mont(a, b, gpu_basic_result, n, p, omega);

    const int gpu_runs = (n <= 16384) ? 3 : 1;
    double gpu_basic_total_time = 0.0;

    for (int run = 0; run < gpu_runs; run++) {
      uint32_t *a_copy = new uint32_t[n];
      uint32_t *b_copy = new uint32_t[n];
      std::copy(a, a + n, a_copy);
      std::copy(b, b + n, b_copy);

      auto gpu_start = high_resolution_clock::now();
      poly_multiply_ntt_gpu_mont(a_copy, b_copy, gpu_basic_result, n, p, omega);
      auto gpu_end = high_resolution_clock::now();

      gpu_basic_total_time +=
          duration_cast<microseconds>(gpu_end - gpu_start).count();

      delete[] a_copy;
      delete[] b_copy;
    }
    result.gpu_basic_time_us = gpu_basic_total_time / gpu_runs;

    // === GPU优化版本测试 ===
    // 预热
    poly_multiply_ntt_gpu_mont_optimized(a, b, gpu_opt_result, n, p, omega);

    double gpu_opt_total_time = 0.0;

    for (int run = 0; run < gpu_runs; run++) {
      uint32_t *a_copy = new uint32_t[n];
      uint32_t *b_copy = new uint32_t[n];
      std::copy(a, a + n, a_copy);
      std::copy(b, b + n, b_copy);

      auto gpu_start = high_resolution_clock::now();
      poly_multiply_ntt_gpu_mont_optimized(a_copy, b_copy, gpu_opt_result, n, p,
                                           omega);
      auto gpu_end = high_resolution_clock::now();

      gpu_opt_total_time +=
          duration_cast<microseconds>(gpu_end - gpu_start).count();

      delete[] a_copy;
      delete[] b_copy;
    }
    result.gpu_opt_time_us = gpu_opt_total_time / gpu_runs;

    // 计算性能指标
    result.speedup_basic = result.cpu_time_us / result.gpu_basic_time_us;
    result.speedup_opt = result.cpu_time_us / result.gpu_opt_time_us;
    result.improvement_ratio =
        result.gpu_basic_time_us / result.gpu_opt_time_us;

    // 验证正确性
    result.correctness_pass =
        verify_results(cpu_result, gpu_basic_result, gpu_opt_result, 2 * n - 1,
                       "n=" + std::to_string(n));

    results.push_back(result);

    // 输出结果
    std::cout << std::setw(8) << result.n << std::setw(12) << std::fixed
              << std::setprecision(1) << result.cpu_time_us << std::setw(12)
              << std::fixed << std::setprecision(1) << result.gpu_basic_time_us
              << std::setw(12) << std::fixed << std::setprecision(1)
              << result.gpu_opt_time_us << std::setw(10) << std::fixed
              << std::setprecision(2) << result.speedup_basic << std::setw(10)
              << std::fixed << std::setprecision(2) << result.speedup_opt
              << std::setw(10) << std::fixed << std::setprecision(2)
              << result.improvement_ratio << std::setw(8)
              << (result.correctness_pass ? "✓" : "✗") << std::endl;

    delete[] a;
    delete[] b;
    delete[] cpu_result;
    delete[] gpu_basic_result;
    delete[] gpu_opt_result;
  }

  // 性能统计分析
  std::cout << std::string(85, '=') << std::endl;
  std::cout << "=== 优化效果统计分析 ===" << std::endl;

  double avg_basic_speedup = 0.0, avg_opt_speedup = 0.0, avg_improvement = 0.0;
  double max_opt_speedup = 0.0;
  uint32_t max_opt_size = 0;
  int valid_count = 0;

  for (const auto &res : results) {
    if (res.correctness_pass) {
      avg_basic_speedup += res.speedup_basic;
      avg_opt_speedup += res.speedup_opt;
      avg_improvement += res.improvement_ratio;
      valid_count++;

      if (res.speedup_opt > max_opt_speedup) {
        max_opt_speedup = res.speedup_opt;
        max_opt_size = res.n;
      }
    }
  }

  if (valid_count > 0) {
    avg_basic_speedup /= valid_count;
    avg_opt_speedup /= valid_count;
    avg_improvement /= valid_count;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "📊 GPU基础版本平均加速比: " << avg_basic_speedup << "x"
              << std::endl;
    std::cout << "🚀 GPU优化版本平均加速比: " << avg_opt_speedup << "x"
              << std::endl;
    std::cout << "📈 优化改进比例: " << avg_improvement << "x" << std::endl;
    std::cout << "🏆 最大优化加速比: " << max_opt_speedup
              << "x (n=" << max_opt_size << ")" << std::endl;

    // 评估优化效果
    std::cout << std::endl << "🎯 优化效果评估: ";
    if (avg_opt_speedup >= 5.0) {
      std::cout << "卓越优化 - 大幅超越理论预期!" << std::endl;
    } else if (avg_opt_speedup >= 3.0) {
      std::cout << "显著优化 - GPU充分发挥并行优势!" << std::endl;
    } else if (avg_opt_speedup >= 2.0) {
      std::cout << "良好优化 - GPU实现了预期加速效果" << std::endl;
    } else if (avg_opt_speedup >= 1.5) {
      std::cout << "一般优化 - 仍有进一步提升空间" << std::endl;
    } else {
      std::cout << "优化不足 - 需要深入分析瓶颈" << std::endl;
    }

    std::cout << std::endl << "相比基础GPU版本的改进: ";
    if (avg_improvement >= 2.0) {
      std::cout << "重大改进 - 优化策略效果显著!" << std::endl;
    } else if (avg_improvement >= 1.5) {
      std::cout << "明显改进 - 优化策略有效" << std::endl;
    } else if (avg_improvement >= 1.2) {
      std::cout << "适度改进 - 优化有一定效果" << std::endl;
    } else {
      std::cout << "改进有限 - 需要更深层次优化" << std::endl;
    }
  }

  // 生成CSV数据用于进一步分析
  std::cout << std::endl << "=== CSV数据 (用于图表生成) ===" << std::endl;
  std::cout << "Size,CPU_Time_us,GPU_Basic_Time_us,GPU_Opt_Time_us,Basic_"
               "Speedup,Opt_Speedup,Improvement_Ratio,Correct"
            << std::endl;
  for (const auto &res : results) {
    std::cout << res.n << "," << res.cpu_time_us << "," << res.gpu_basic_time_us
              << "," << res.gpu_opt_time_us << "," << res.speedup_basic << ","
              << res.speedup_opt << "," << res.improvement_ratio << ","
              << (res.correctness_pass ? "true" : "false") << std::endl;
  }
}

int main() {
  // 检查CUDA设备
  int deviceCount;
  cudaError_t error = cudaGetDeviceCount(&deviceCount);

  if (error != cudaSuccess || deviceCount == 0) {
    std::cerr << "❌ 错误: 未检测到CUDA设备!" << std::endl;
    std::cerr << "请检查GPU驱动和CUDA安装" << std::endl;
    return 1;
  }

  // 显示GPU信息
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  std::cout << "🔧 GPU设备: " << prop.name << std::endl;
  std::cout << "🔧 计算能力: " << prop.major << "." << prop.minor << std::endl;
  std::cout << "🔧 SM数量: " << prop.multiProcessorCount << std::endl;
  std::cout << "🔧 全局内存: " << (prop.totalGlobalMem >> 20) << " MB"
            << std::endl;
  std::cout << "🔧 最大线程/块: " << prop.maxThreadsPerBlock << std::endl;
  std::cout << std::endl;

  try {
    run_comprehensive_benchmark();
  } catch (const std::exception &e) {
    std::cerr << "❌ 测试错误: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}