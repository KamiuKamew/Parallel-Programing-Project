#include "src/include/CUDA/ntt.h"
#include "src/include/ntt.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// CUDA精确计时器
class CudaTimer {
private:
  cudaEvent_t start_event, stop_event;

public:
  CudaTimer() {
    cudaEventCreate(&start_event);
    cudaEventCreate(&stop_event);
  }

  ~CudaTimer() {
    cudaEventDestroy(start_event);
    cudaEventDestroy(stop_event);
  }

  void start() {
    cudaDeviceSynchronize(); // 确保之前的操作完成
    cudaEventRecord(start_event);
  }

  double stop() {
    cudaEventRecord(stop_event);
    cudaEventSynchronize(stop_event);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start_event, stop_event);
    return milliseconds * 1000.0; // 转换为微秒
  }
};

// CPU精确计时器
class CpuTimer {
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

// 单个算法性能测试结果
template <typename T> struct AlgorithmResults {
  double cpu_time;
  double gpu_time;
  double speedup;
  bool correctness;
  std::string algorithm_name;
};

// GPU热身运行
template <typename T> void gpu_warmup(T n, T p, T omega) {
  T *a = new T[n];
  T *b = new T[n];
  T *ab = new T[2 * n - 1];

  generate_test_case(a, b, n, 0);

  // 各算法热身一次
  poly_multiply_ntt_gpu_naive(a, b, ab, n, p, omega);
  cudaDeviceSynchronize();

  poly_multiply_ntt_gpu_mont(a, b, ab, n, p, omega);
  cudaDeviceSynchronize();

  poly_multiply_ntt_gpu_barrett(a, b, ab, n, p, omega);
  cudaDeviceSynchronize();

  delete[] a;
  delete[] b;
  delete[] ab;

  std::cout << "GPU warmup completed for n=" << n << std::endl;
}

// 执行单个模乘算法测试（修复版）
template <typename T>
AlgorithmResults<T> test_single_algorithm_fixed(T n, T p, T omega,
                                                const std::string &algo_name,
                                                int mode) {
  AlgorithmResults<T> result;
  result.algorithm_name = algo_name;

  // 分配内存
  T *a = new T[n];
  T *b = new T[n];
  T *cpu_result = new T[2 * n - 1];
  T *gpu_result = new T[2 * n - 1];

  // 生成固定测试数据
  generate_test_case(a, b, n, 0);

  CpuTimer cpu_timer;
  CudaTimer gpu_timer;
  const int num_runs = (n <= 1024) ? 5 : 3;

  // 测试CPU版本（只有Montgomery可用）
  double cpu_total_time = 0.0;
  if (algo_name == "Montgomery") {
    for (int run = 0; run < num_runs; run++) {
      // 每次使用相同的测试数据以保证公平比较
      T *a_copy = new T[n];
      T *b_copy = new T[n];
      std::copy(a, a + n, a_copy);
      std::copy(b, b + n, b_copy);

      cpu_timer.start();
      poly_multiply_ntt(a_copy, b_copy, cpu_result, n, p, omega);
      double cpu_time = cpu_timer.stop();
      cpu_total_time += cpu_time;

      delete[] a_copy;
      delete[] b_copy;
    }
    result.cpu_time = cpu_total_time / num_runs;
  } else {
    result.cpu_time = -1; // 标记不可用
  }

  // 测试GPU版本（使用CUDA计时器进行精确计时）
  double gpu_total_time = 0.0;
  for (int run = 0; run < num_runs; run++) {
    // 每次使用相同的测试数据
    T *a_copy = new T[n];
    T *b_copy = new T[n];
    std::copy(a, a + n, a_copy);
    std::copy(b, b + n, b_copy);

    gpu_timer.start();
    switch (mode) {
    case 0: // Naive
      poly_multiply_ntt_gpu_naive(a_copy, b_copy, gpu_result, n, p, omega);
      break;
    case 1: // Montgomery
      poly_multiply_ntt_gpu_mont(a_copy, b_copy, gpu_result, n, p, omega);
      break;
    case 2: // Barrett
      poly_multiply_ntt_gpu_barrett(a_copy, b_copy, gpu_result, n, p, omega);
      break;
    }
    double gpu_time = gpu_timer.stop();
    gpu_total_time += gpu_time;

    delete[] a_copy;
    delete[] b_copy;
  }
  result.gpu_time = gpu_total_time / num_runs;

  // 计算加速比
  if (result.cpu_time > 0) {
    result.speedup = result.cpu_time / result.gpu_time;
  } else {
    result.speedup = -1; // 无法计算
  }

  // 验证正确性（使用固定数据）
  T *a_test = new T[n];
  T *b_test = new T[n];
  std::copy(a, a + n, a_test);
  std::copy(b, b + n, b_test);

  poly_multiply_ntt(a_test, b_test, cpu_result, n, p, omega); // CPU基线

  std::copy(a, a + n, a_test);
  std::copy(b, b + n, b_test);

  switch (mode) {
  case 0:
    poly_multiply_ntt_gpu_naive(a_test, b_test, gpu_result, n, p, omega);
    break;
  case 1:
    poly_multiply_ntt_gpu_mont(a_test, b_test, gpu_result, n, p, omega);
    break;
  case 2:
    poly_multiply_ntt_gpu_barrett(a_test, b_test, gpu_result, n, p, omega);
    break;
  }

  result.correctness = verify_results(cpu_result, gpu_result, 2 * n - 1);

  // 清理内存
  delete[] a;
  delete[] b;
  delete[] cpu_result;
  delete[] gpu_result;
  delete[] a_test;
  delete[] b_test;

  return result;
}

// 执行完整的性能对比测试（修复版）
template <typename T>
void run_comprehensive_test_fixed(T n, u32 modulus, u32 omega,
                                  const std::string &modulus_name) {
  std::cout << "\n=== Fixed Performance Test (n=" << n << ", p=" << modulus
            << ") ===" << std::endl;
  std::cout << "Modulus: " << modulus_name << std::endl;

  // GPU热身
  gpu_warmup<T>(n, modulus, omega);

  std::vector<AlgorithmResults<T>> results;

  // 测试三种模乘算法
  std::cout << "Testing Naive algorithm..." << std::endl;
  results.push_back(
      test_single_algorithm_fixed<T>(n, modulus, omega, "Naive", 0));

  std::cout << "Testing Montgomery algorithm..." << std::endl;
  results.push_back(
      test_single_algorithm_fixed<T>(n, modulus, omega, "Montgomery", 1));

  std::cout << "Testing Barrett algorithm..." << std::endl;
  results.push_back(
      test_single_algorithm_fixed<T>(n, modulus, omega, "Barrett", 2));

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

  std::cout << "Montgomery vs Naive: " << std::fixed << std::setprecision(2)
            << (naive_gpu / mont_gpu) << "x speedup" << std::endl;
  std::cout << "Barrett vs Naive: " << (naive_gpu / barrett_gpu) << "x speedup"
            << std::endl;
  std::cout << "Montgomery vs Barrett: " << (barrett_gpu / mont_gpu) << "x"
            << ((mont_gpu < barrett_gpu) ? " speedup" : " slower") << std::endl;

  // 输出到CSV
  std::ofstream csv_file("fixed_modmul_results.csv", std::ios::app);
  if (csv_file.tellp() == 0) {
    csv_file << "Algorithm,Platform,N,Modulus,Time_us,Correctness,Speedup_vs_"
                "CPU_Naive,Speedup_vs_GPU_Naive\n";
  }

  for (const auto &result : results) {
    if (result.cpu_time > 0) {
      csv_file << "CPU " << result.algorithm_name << ",CPU," << n << ","
               << modulus << "," << std::fixed << std::setprecision(3)
               << result.cpu_time << ","
               << (result.correctness ? "PASS" : "FAIL") << ",1.000,0.000\n";
    }

    csv_file << "GPU " << result.algorithm_name << ",GPU," << n << ","
             << modulus << "," << std::fixed << std::setprecision(3)
             << result.gpu_time << "," << (result.correctness ? "PASS" : "FAIL")
             << ",";

    if (result.cpu_time > 0) {
      csv_file << result.speedup;
    } else {
      csv_file << "0.000";
    }

    csv_file << "," << (naive_gpu / result.gpu_time) << "\n";
  }
  csv_file.close();
}

int main() {
  std::cout << "=== Fixed Comprehensive Modular Multiplication Performance "
               "Analysis ==="
            << std::endl;
  std::cout << "Testing with CUDA event timing and proper memory management"
            << std::endl;

  // 定义测试配置
  struct TestConfig {
    u32 modulus;
    u32 omega;
    std::string name;
  };

  std::vector<TestConfig> test_configs = {
      {7340033u, 3u, "Small Modulus (7.3M)"},
      {104857601u, 3u, "Medium Modulus (104M)"},
      {469762049u, 3u, "Large Modulus (469M)"},
      {998244353u, 3u, "Large Modulus (998M)"},
  };

  // 测试不同规模
  std::vector<u32> test_sizes = {4, 1024, 16384, 65536, 131072};

  // 清空之前的结果文件
  std::ofstream clear_file("fixed_modmul_results.csv", std::ios::trunc);
  clear_file.close();

  for (const auto &config : test_configs) {
    for (u32 n : test_sizes) {
      // 对于大规模测试，跳过小模数以节省时间
      if (n > 16384 && config.modulus < 100000000u) {
        continue;
      }

      run_comprehensive_test_fixed<u32>(n, config.modulus, config.omega,
                                        config.name);
    }
  }

  std::cout << "\n=== Fixed Performance Analysis Complete ===" << std::endl;
  std::cout << "Results saved to fixed_modmul_results.csv" << std::endl;

  return 0;
}