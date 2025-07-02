#include "src/include/CUDA/ntt_gpu_optimized.cu"
#include "src/include/ntt.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std::chrono;

// 性能测试结果结构
struct BenchmarkResult {
  uint32_t n;
  double cpu_time_us;
  double gpu_basic_time_us;
  double gpu_optimized_time_us;
  double speedup_basic;
  double speedup_optimized;
  double improvement_ratio;
};

// 测试数据生成
template <typename T> void generate_test_data(T *a, T *b, T n, T p) {
  for (T i = 0; i < n; i++) {
    a[i] = (i + 1) % p;
    b[i] = (i * 2 + 1) % p;
  }
}

// 结果验证
template <typename T> bool verify_results(T *cpu_result, T *gpu_result, T n) {
  for (T i = 0; i < 2 * n - 1; i++) {
    if (cpu_result[i] != gpu_result[i]) {
      std::cout << "Verification failed at index " << i
                << ": CPU=" << cpu_result[i] << ", GPU=" << gpu_result[i]
                << std::endl;
      return false;
    }
  }
  return true;
}

// 单次性能测试
template <typename T>
BenchmarkResult benchmark_single_case(T n, T p, T omega, int warmup_runs = 3,
                                      int test_runs = 10) {
  BenchmarkResult result;
  result.n = n;

  // 分配测试数据
  T *a = new T[n];
  T *b = new T[n];
  T *cpu_result = new T[2 * n - 1];
  T *gpu_basic_result = new T[2 * n - 1];
  T *gpu_optimized_result = new T[2 * n - 1];

  generate_test_data(a, b, n, p);

  std::cout << "Testing n=" << n << "..." << std::endl;

  // 预热GPU
  for (int i = 0; i < warmup_runs; i++) {
    poly_multiply_ntt_gpu_mont(a, b, gpu_basic_result, n, p, omega);
    poly_multiply_ntt_gpu_mont_optimized(a, b, gpu_optimized_result, n, p,
                                         omega);
  }

  // CPU基准测试
  auto cpu_start = high_resolution_clock::now();
  for (int i = 0; i < test_runs; i++) {
    poly_multiply_ntt(a, b, cpu_result, n, p, omega);
  }
  auto cpu_end = high_resolution_clock::now();
  result.cpu_time_us =
      duration_cast<microseconds>(cpu_end - cpu_start).count() / test_runs;

  // GPU基础版本测试
  auto gpu_basic_start = high_resolution_clock::now();
  for (int i = 0; i < test_runs; i++) {
    poly_multiply_ntt_gpu_mont(a, b, gpu_basic_result, n, p, omega);
  }
  auto gpu_basic_end = high_resolution_clock::now();
  result.gpu_basic_time_us =
      duration_cast<microseconds>(gpu_basic_end - gpu_basic_start).count() /
      test_runs;

  // GPU优化版本测试
  auto gpu_opt_start = high_resolution_clock::now();
  for (int i = 0; i < test_runs; i++) {
    poly_multiply_ntt_gpu_mont_optimized(a, b, gpu_optimized_result, n, p,
                                         omega);
  }
  auto gpu_opt_end = high_resolution_clock::now();
  result.gpu_optimized_time_us =
      duration_cast<microseconds>(gpu_opt_end - gpu_opt_start).count() /
      test_runs;

  // 验证正确性
  bool cpu_vs_basic = verify_results(cpu_result, gpu_basic_result, n);
  bool cpu_vs_optimized = verify_results(cpu_result, gpu_optimized_result, n);

  if (!cpu_vs_basic) {
    std::cout << "ERROR: GPU basic version produces incorrect results!"
              << std::endl;
  }
  if (!cpu_vs_optimized) {
    std::cout << "ERROR: GPU optimized version produces incorrect results!"
              << std::endl;
  }

  // 计算性能指标
  result.speedup_basic = result.cpu_time_us / result.gpu_basic_time_us;
  result.speedup_optimized = result.cpu_time_us / result.gpu_optimized_time_us;
  result.improvement_ratio =
      result.gpu_basic_time_us / result.gpu_optimized_time_us;

  delete[] a;
  delete[] b;
  delete[] cpu_result;
  delete[] gpu_basic_result;
  delete[] gpu_optimized_result;

  return result;
}

// 完整性能测试套件
void run_comprehensive_benchmark() {
  const u32 p = 998244353; // 常用NTT模数
  const u32 omega = 3;     // 原根

  // 测试规模范围：从4到1M
  std::vector<u32> test_sizes = {4,     16,    64,     256,    1024,   4096,
                                 16384, 65536, 131072, 262144, 524288, 1048576};

  std::vector<BenchmarkResult> results;

  std::cout << "=== GPU NTT优化效果综合测试 ===" << std::endl;
  std::cout << std::setw(10) << "Size" << std::setw(12) << "CPU(μs)"
            << std::setw(12) << "GPU基础(μs)" << std::setw(12) << "GPU优化(μs)"
            << std::setw(12) << "基础加速比" << std::setw(12) << "优化加速比"
            << std::setw(12) << "优化提升" << std::endl;
  std::cout << std::string(88, '-') << std::endl;

  for (u32 n : test_sizes) {
    BenchmarkResult result = benchmark_single_case<u32>(n, p, omega);
    results.push_back(result);

    std::cout << std::setw(10) << result.n << std::setw(12) << std::fixed
              << std::setprecision(1) << result.cpu_time_us << std::setw(12)
              << std::fixed << std::setprecision(1) << result.gpu_basic_time_us
              << std::setw(12) << std::fixed << std::setprecision(1)
              << result.gpu_optimized_time_us << std::setw(12) << std::fixed
              << std::setprecision(3) << result.speedup_basic << std::setw(12)
              << std::fixed << std::setprecision(3) << result.speedup_optimized
              << std::setw(12) << std::fixed << std::setprecision(3)
              << result.improvement_ratio << std::endl;
  }

  // 保存结果到CSV文件
  std::ofstream csv_file("gpu_optimization_benchmark.csv");
  csv_file << "Size,CPU_Time_us,GPU_Basic_Time_us,GPU_Optimized_Time_us,Basic_"
              "Speedup,Optimized_Speedup,Improvement_Ratio\n";

  for (const auto &result : results) {
    csv_file << result.n << "," << result.cpu_time_us << ","
             << result.gpu_basic_time_us << "," << result.gpu_optimized_time_us
             << "," << result.speedup_basic << "," << result.speedup_optimized
             << "," << result.improvement_ratio << "\n";
  }
  csv_file.close();

  std::cout << "\n结果已保存到 gpu_optimization_benchmark.csv" << std::endl;

  // 输出总结分析
  std::cout << "\n=== 优化效果分析 ===" << std::endl;

  // 找到最大加速比
  auto max_speedup_basic =
      std::max_element(results.begin(), results.end(),
                       [](const BenchmarkResult &a, const BenchmarkResult &b) {
                         return a.speedup_basic < b.speedup_basic;
                       });

  auto max_speedup_optimized =
      std::max_element(results.begin(), results.end(),
                       [](const BenchmarkResult &a, const BenchmarkResult &b) {
                         return a.speedup_optimized < b.speedup_optimized;
                       });

  auto max_improvement =
      std::max_element(results.begin(), results.end(),
                       [](const BenchmarkResult &a, const BenchmarkResult &b) {
                         return a.improvement_ratio < b.improvement_ratio;
                       });

  std::cout << "最佳基础GPU加速比: " << max_speedup_basic->speedup_basic
            << "x (n=" << max_speedup_basic->n << ")" << std::endl;
  std::cout << "最佳优化GPU加速比: " << max_speedup_optimized->speedup_optimized
            << "x (n=" << max_speedup_optimized->n << ")" << std::endl;
  std::cout << "最大优化提升倍数: " << max_improvement->improvement_ratio
            << "x (n=" << max_improvement->n << ")" << std::endl;

  // 计算大规模问题平均提升
  double avg_improvement_large = 0;
  int large_count = 0;
  for (const auto &result : results) {
    if (result.n >= 65536) { // 大规模问题
      avg_improvement_large += result.improvement_ratio;
      large_count++;
    }
  }
  if (large_count > 0) {
    avg_improvement_large /= large_count;
    std::cout << "大规模问题平均优化提升: " << avg_improvement_large << "x"
              << std::endl;
  }
}

// 详细性能分析：测试不同优化策略的独立效果
void analyze_optimization_strategies() {
  const u32 n = 131072; // 中等规模测试
  const u32 p = 998244353;
  const u32 omega = 3;

  std::cout << "\n=== 优化策略独立效果分析 (n=" << n << ") ===" << std::endl;

  u32 *a = new u32[n];
  u32 *b = new u32[n];
  u32 *result = new u32[2 * n - 1];

  generate_test_data(a, b, n, p);

  // 测试基础版本
  auto start = high_resolution_clock::now();
  poly_multiply_ntt_gpu_mont(a, b, result, n, p, omega);
  auto end = high_resolution_clock::now();
  double basic_time = duration_cast<microseconds>(end - start).count();

  // 测试优化版本
  start = high_resolution_clock::now();
  poly_multiply_ntt_gpu_mont_optimized(a, b, result, n, p, omega);
  end = high_resolution_clock::now();
  double optimized_time = duration_cast<microseconds>(end - start).count();

  std::cout << "基础版本耗时: " << basic_time << " μs" << std::endl;
  std::cout << "优化版本耗时: " << optimized_time << " μs" << std::endl;
  std::cout << "性能提升: " << (basic_time / optimized_time) << "x"
            << std::endl;

  delete[] a;
  delete[] b;
  delete[] result;
}

int main() {
  std::cout << "GPU NTT优化效果测试程序" << std::endl;
  std::cout << "========================" << std::endl;

  // 获取GPU信息
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  std::cout << "GPU: " << prop.name << std::endl;
  std::cout << "SM数量: " << prop.multiProcessorCount << std::endl;
  std::cout << "全局内存: " << prop.totalGlobalMem / (1024 * 1024) << " MB"
            << std::endl;
  std::cout << "共享内存/块: " << prop.sharedMemPerBlock / 1024 << " KB"
            << std::endl;
  std::cout << std::endl;

  try {
    // 运行综合性能测试
    run_comprehensive_benchmark();

    // 运行详细策略分析
    analyze_optimization_strategies();

  } catch (const std::exception &e) {
    std::cerr << "测试过程中发生错误: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}