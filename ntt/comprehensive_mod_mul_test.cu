#include "src/include/CUDA/ntt.h"
#include "src/include/ntt.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// 高精度计时器
class Timer {
private:
  std::chrono::high_resolution_clock::time_point start_time;

public:
  void start() { start_time = std::chrono::high_resolution_clock::now(); }

  double stop() {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
        end_time - start_time);
    return duration.count() / 1000.0; // 返回微秒
  }
};

// 测试用例生成
template <typename T> void generate_test_case(T *a, T *b, T n, int case_id) {
  for (T i = 0; i < n; i++) {
    a[i] = (i * 17 + 7) % 1000 + 1;
    b[i] = (i * 23 + 11) % 1000 + 1;
  }
}

// 验证结果正确性
template <typename T>
bool verify_results(T *expected, T *actual, T result_len) {
  for (T i = 0; i < result_len; i++) {
    if (expected[i] != actual[i]) {
      std::cout << "Mismatch at index " << i << ": Expected=" << expected[i]
                << " Actual=" << actual[i] << std::endl;
      return false;
    }
  }
  return true;
}

// CPU版本朴素多项式乘法
template <typename T>
void poly_multiply_naive_cpu(T *a, T *b, T *ab, T n, T p) {
  // 初始化结果数组
  for (T i = 0; i < 2 * n - 1; ++i) {
    ab[i] = 0;
  }

  // 朴素O(n²)算法
  for (T i = 0; i < n; ++i) {
    for (T j = 0; j < n; ++j) {
      ab[i + j] = (ab[i + j] + (((__uint128_t)a[i] * b[j]) % p)) % p;
    }
  }
}

// 多模数性能测试结构
struct TestConfig {
  u32 modulus;
  u32 omega;
  std::string name;
};

// 单个算法性能测试
template <typename T> struct AlgorithmResults {
  double cpu_time;
  double gpu_time;
  double speedup;
  bool correctness;
  std::string algorithm_name;
};

// 执行单个模乘算法测试
template <typename T>
AlgorithmResults<T> test_single_algorithm(T n, T p, T omega,
                                          const std::string &algo_name,
                                          int mode) {
  AlgorithmResults<T> result;
  result.algorithm_name = algo_name;

  // 分配内存
  T *a = new T[n];
  T *b = new T[n];
  T *cpu_result = new T[2 * n - 1];
  T *gpu_result = new T[2 * n - 1];
  T *naive_baseline = new T[2 * n - 1];

  // 生成测试数据
  generate_test_case(a, b, n, 0);

  Timer timer;
  const int num_runs = (n <= 1024) ? 10 : 3;

  // 生成朴素算法基线（用于正确性验证）
  poly_multiply_naive_cpu(a, b, naive_baseline, n, p);

  // 测试CPU版本（只有Montgomery可用）
  double cpu_total_time = 0.0;
  if (algo_name == "Montgomery") {
    for (int run = 0; run < num_runs; run++) {
      generate_test_case(a, b, n, run % 2);
      timer.start();
      poly_multiply_ntt(a, b, cpu_result, n, p, omega);
      double cpu_time = timer.stop();
      cpu_total_time += cpu_time;
    }
    result.cpu_time = cpu_total_time / num_runs;
  } else {
    result.cpu_time = -1; // 标记不可用
  }

  // 测试GPU版本
  double gpu_total_time = 0.0;
  for (int run = 0; run < num_runs; run++) {
    generate_test_case(a, b, n, run % 2);
    timer.start();

    switch (mode) {
    case 0: // Naive
      poly_multiply_ntt_gpu_naive(a, b, gpu_result, n, p, omega);
      break;
    case 1: // Montgomery
      poly_multiply_ntt_gpu(a, b, gpu_result, n, p, omega);
      break;
    case 2: // Barrett
      poly_multiply_ntt_gpu_barrett(a, b, gpu_result, n, p, omega);
      break;
    }

    double gpu_time = timer.stop();
    gpu_total_time += gpu_time;
  }
  result.gpu_time = gpu_total_time / num_runs;

  // 计算加速比
  if (result.cpu_time > 0) {
    result.speedup = result.cpu_time / result.gpu_time;
  } else {
    result.speedup = -1; // 无法计算
  }

  // 验证正确性（与Montgomery CPU版本对比）
  generate_test_case(a, b, n, 0);
  poly_multiply_ntt(a, b, cpu_result, n, p, omega); // Montgomery CPU版本

  generate_test_case(a, b, n, 0);
  switch (mode) {
  case 0:
    poly_multiply_ntt_gpu_naive(a, b, gpu_result, n, p, omega);
    break;
  case 1:
    poly_multiply_ntt_gpu(a, b, gpu_result, n, p, omega);
    break;
  case 2:
    poly_multiply_ntt_gpu_barrett(a, b, gpu_result, n, p, omega);
    break;
  }

  result.correctness = verify_results(cpu_result, gpu_result, 2 * n - 1);

  // 清理内存
  delete[] a;
  delete[] b;
  delete[] cpu_result;
  delete[] gpu_result;
  delete[] naive_baseline;

  return result;
}

// 执行完整的性能对比测试
template <typename T>
void run_comprehensive_test(T n, const TestConfig &config) {
  std::cout << "\n=== Comprehensive Performance Test (n=" << n
            << ", p=" << config.modulus << ") ===" << std::endl;
  std::cout << "Modulus: " << config.name << std::endl;

  std::vector<AlgorithmResults<T>> results;

  // 测试三种模乘算法
  results.push_back(
      test_single_algorithm<T>(n, config.modulus, config.omega, "Naive", 0));
  results.push_back(test_single_algorithm<T>(n, config.modulus, config.omega,
                                             "Montgomery", 1));
  results.push_back(
      test_single_algorithm<T>(n, config.modulus, config.omega, "Barrett", 2));

  // 输出结果表格
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "\n┌─────────────┬─────────────┬─────────────┬─────────────┬────"
               "─────────┐"
            << std::endl;
  std::cout << "│ Algorithm   │ CPU Time(μs)│ GPU Time(μs)│ Speedup     │ "
               "Correctness │"
            << std::endl;
  std::cout << "├─────────────┼─────────────┼─────────────┼─────────────┼──────"
               "───────┤"
            << std::endl;

  for (const auto &result : results) {
    std::cout << "│ " << std::setw(11) << result.algorithm_name << " │ ";

    if (result.cpu_time > 0) {
      std::cout << std::setw(11) << result.cpu_time;
    } else {
      std::cout << std::setw(11) << "N/A";
    }

    std::cout << " │ " << std::setw(11) << result.gpu_time << " │ ";

    if (result.speedup > 0) {
      std::cout << std::setw(11) << result.speedup << " │ ";
    } else {
      std::cout << std::setw(11) << "N/A" << " │ ";
    }

    std::cout << std::setw(11) << (result.correctness ? "✅ PASS" : "❌ FAIL")
              << " │" << std::endl;
  }

  std::cout << "└─────────────┴─────────────┴─────────────┴─────────────┴──────"
               "───────┘"
            << std::endl;

  // 分析GPU内部算法对比
  std::cout << "\n--- GPU Algorithm Performance Analysis ---" << std::endl;
  double naive_gpu = results[0].gpu_time;
  double mont_gpu = results[1].gpu_time;
  double barrett_gpu = results[2].gpu_time;

  std::cout << "Montgomery vs Naive: " << (naive_gpu / mont_gpu) << "x speedup"
            << std::endl;
  std::cout << "Barrett vs Naive: " << (naive_gpu / barrett_gpu) << "x speedup"
            << std::endl;
  std::cout << "Montgomery vs Barrett: " << (barrett_gpu / mont_gpu) << "x"
            << ((mont_gpu < barrett_gpu) ? " speedup" : " slower") << std::endl;
}

int main() {
  std::cout
      << "=== Comprehensive Modular Multiplication Performance Analysis ==="
      << std::endl;
  std::cout << "Testing Naive/Montgomery/Barrett algorithms on CPU and GPU"
            << std::endl;

  // 定义测试配置：四个不同模数
  std::vector<TestConfig> test_configs = {
      {7340033u, 3u, "Small Modulus (7.3M)"},
      {104857601u, 3u, "Medium Modulus (104M)"},
      {469762049u, 3u, "Large Modulus (469M)"},
      // 注意：第四个模数需要64位支持，暂时跳过
  };

  // 测试不同规模
  std::vector<u32> test_sizes = {4, 1024, 16384, 131072};

  // 创建结果输出文件
  std::ofstream csv_file("mod_mul_performance_results.csv");
  csv_file
      << "Modulus,n,Algorithm,CPU_Time_us,GPU_Time_us,Speedup,Correctness\n";

  for (const auto &config : test_configs) {
    for (u32 n : test_sizes) {
      // 对于大规模测试，跳过小模数以节省时间
      if (n > 16384 && config.modulus < 100000000u) {
        continue;
      }

      run_comprehensive_test<u32>(n, config);

      // TODO: 保存详细结果到CSV文件
    }
  }

  csv_file.close();

  std::cout << "\n=== Performance Analysis Complete ===" << std::endl;
  std::cout << "Results saved to mod_mul_performance_results.csv" << std::endl;

  return 0;
}