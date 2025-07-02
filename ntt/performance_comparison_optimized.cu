#include "src/include/CUDA/ntt.h"
#include "src/include/CUDA/ntt_gpu_optimized.cu"
#include "src/include/ntt.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

// 高精度CUDA计时器
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
    cudaDeviceSynchronize();
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

// CPU计时器
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

// 性能测试结果结构
struct PerformanceResult {
  std::string algorithm_name;
  std::string implementation_type; // "CPU", "GPU_Original", "GPU_Optimized"
  double execution_time_us;
  double speedup_vs_cpu;
  double speedup_vs_original_gpu;
  bool correctness_verified;
};

// 测试数据生成
template <typename T> void generate_test_case(T *a, T *b, T n) {
  for (T i = 0; i < n; i++) {
    a[i] = (i * 17 + 7) % 1000 + 1;
    b[i] = (i * 23 + 11) % 1000 + 1;
  }
}

// 结果验证
template <typename T>
bool verify_results(T *reference, T *test, T result_len,
                    const std::string &test_name) {
  for (T i = 0; i < result_len; i++) {
    if (reference[i] != test[i]) {
      std::cout << "❌ " << test_name << " 验证失败: index=" << i
                << " reference=" << reference[i] << " test=" << test[i]
                << std::endl;
      return false;
    }
  }
  std::cout << "✅ " << test_name << " 验证通过" << std::endl;
  return true;
}

// GPU热身
template <typename T> void gpu_warmup(T n, T p, T omega) {
  T *a = new T[n];
  T *b = new T[n];
  T *ab = new T[2 * n - 1];

  generate_test_case(a, b, n);

  // 原始版本热身
  poly_multiply_ntt_gpu_mont(a, b, ab, n, p, omega);
  cudaDeviceSynchronize();

  // 优化版本热身
  poly_multiply_ntt_gpu_mont_optimized(a, b, ab, n, p, omega);
  cudaDeviceSynchronize();

  delete[] a;
  delete[] b;
  delete[] ab;
}

// 单个算法的完整性能测试
template <typename T>
std::vector<PerformanceResult>
test_algorithm_comprehensive(T n, T p, T omega, const std::string &algo_name) {
  std::vector<PerformanceResult> results;

  std::cout << "\n=== 测试 " << algo_name << " 算法 (n=" << n << ", p=" << p
            << ") ===" << std::endl;

  // 分配内存
  T *a = new T[n];
  T *b = new T[n];
  T *cpu_result = new T[2 * n - 1];
  T *gpu_orig_result = new T[2 * n - 1];
  T *gpu_opt_result = new T[2 * n - 1];

  // 生成固定测试数据
  generate_test_case(a, b, n);

  CpuTimer cpu_timer;
  CudaTimer gpu_timer;
  const int num_runs = (n <= 1024) ? 5 : 3;

  // GPU热身
  gpu_warmup(n, p, omega);

  // ========== CPU版本测试 ==========
  double cpu_total_time = 0.0;
  if (algo_name == "Montgomery") { // 只有Montgomery有CPU版本
    for (int run = 0; run < num_runs; run++) {
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

    PerformanceResult cpu_result_entry;
    cpu_result_entry.algorithm_name = algo_name;
    cpu_result_entry.implementation_type = "CPU";
    cpu_result_entry.execution_time_us = cpu_total_time / num_runs;
    cpu_result_entry.speedup_vs_cpu = 1.0;
    cpu_result_entry.speedup_vs_original_gpu = -1; // 不适用
    cpu_result_entry.correctness_verified = true;
    results.push_back(cpu_result_entry);

    std::cout << "CPU " << algo_name << ": " << std::fixed
              << std::setprecision(1) << cpu_result_entry.execution_time_us
              << " μs" << std::endl;
  }

  // ========== GPU原始版本测试 ==========
  double gpu_orig_total_time = 0.0;
  for (int run = 0; run < num_runs; run++) {
    T *a_copy = new T[n];
    T *b_copy = new T[n];
    std::copy(a, a + n, a_copy);
    std::copy(b, b + n, b_copy);

    gpu_timer.start();
    if (algo_name == "Naive") {
      poly_multiply_ntt_gpu_naive(a_copy, b_copy, gpu_orig_result, n, p, omega);
    } else if (algo_name == "Montgomery") {
      poly_multiply_ntt_gpu_mont(a_copy, b_copy, gpu_orig_result, n, p, omega);
    } else if (algo_name == "Barrett") {
      poly_multiply_ntt_gpu_barrett(a_copy, b_copy, gpu_orig_result, n, p,
                                    omega);
    }
    double gpu_time = gpu_timer.stop();
    gpu_orig_total_time += gpu_time;

    delete[] a_copy;
    delete[] b_copy;
  }

  PerformanceResult gpu_orig_result_entry;
  gpu_orig_result_entry.algorithm_name = algo_name;
  gpu_orig_result_entry.implementation_type = "GPU_Original";
  gpu_orig_result_entry.execution_time_us = gpu_orig_total_time / num_runs;
  gpu_orig_result_entry.speedup_vs_cpu =
      (results.empty()) ? -1
                        : results[0].execution_time_us /
                              gpu_orig_result_entry.execution_time_us;
  gpu_orig_result_entry.speedup_vs_original_gpu = 1.0;

  // 验证正确性（与CPU对比，如果CPU版本存在）
  if (!results.empty()) {
    gpu_orig_result_entry.correctness_verified =
        verify_results(cpu_result, gpu_orig_result, 2 * n - 1, "GPU原始版本");
  } else {
    gpu_orig_result_entry.correctness_verified = true; // 无法验证，假设正确
  }

  results.push_back(gpu_orig_result_entry);

  std::cout << "GPU原始 " << algo_name << ": " << std::fixed
            << std::setprecision(1) << gpu_orig_result_entry.execution_time_us
            << " μs";
  if (gpu_orig_result_entry.speedup_vs_cpu > 0) {
    std::cout << " (加速比: " << std::setprecision(2)
              << gpu_orig_result_entry.speedup_vs_cpu << "x)";
  }
  std::cout << std::endl;

  // ========== GPU优化版本测试 ==========
  double gpu_opt_total_time = 0.0;
  for (int run = 0; run < num_runs; run++) {
    T *a_copy = new T[n];
    T *b_copy = new T[n];
    std::copy(a, a + n, a_copy);
    std::copy(b, b + n, b_copy);

    gpu_timer.start();
    if (algo_name == "Naive") {
      poly_multiply_ntt_gpu_naive_optimized(a_copy, b_copy, gpu_opt_result, n,
                                            p, omega);
    } else if (algo_name == "Montgomery") {
      poly_multiply_ntt_gpu_mont_optimized(a_copy, b_copy, gpu_opt_result, n, p,
                                           omega);
    } else if (algo_name == "Barrett") {
      poly_multiply_ntt_gpu_barrett_optimized(a_copy, b_copy, gpu_opt_result, n,
                                              p, omega);
    }
    double gpu_time = gpu_timer.stop();
    gpu_opt_total_time += gpu_time;

    delete[] a_copy;
    delete[] b_copy;
  }

  PerformanceResult gpu_opt_result_entry;
  gpu_opt_result_entry.algorithm_name = algo_name;
  gpu_opt_result_entry.implementation_type = "GPU_Optimized";
  gpu_opt_result_entry.execution_time_us = gpu_opt_total_time / num_runs;
  gpu_opt_result_entry.speedup_vs_cpu =
      (results.empty()) ? -1
                        : results[0].execution_time_us /
                              gpu_opt_result_entry.execution_time_us;
  gpu_opt_result_entry.speedup_vs_original_gpu =
      gpu_orig_result_entry.execution_time_us /
      gpu_opt_result_entry.execution_time_us;

  // 验证正确性
  if (!results.empty()) {
    gpu_opt_result_entry.correctness_verified =
        verify_results(cpu_result, gpu_opt_result, 2 * n - 1, "GPU优化版本");
  } else {
    // 与原始GPU版本对比
    gpu_opt_result_entry.correctness_verified = verify_results(
        gpu_orig_result, gpu_opt_result, 2 * n - 1, "GPU优化版本");
  }

  results.push_back(gpu_opt_result_entry);

  std::cout << "GPU优化 " << algo_name << ": " << std::fixed
            << std::setprecision(1) << gpu_opt_result_entry.execution_time_us
            << " μs";
  if (gpu_opt_result_entry.speedup_vs_cpu > 0) {
    std::cout << " (vs CPU: " << std::setprecision(2)
              << gpu_opt_result_entry.speedup_vs_cpu << "x)";
  }
  std::cout << " (vs 原始GPU: " << std::setprecision(2)
            << gpu_opt_result_entry.speedup_vs_original_gpu << "x)"
            << std::endl;

  // 清理内存
  delete[] a;
  delete[] b;
  delete[] cpu_result;
  delete[] gpu_orig_result;
  delete[] gpu_opt_result;

  return results;
}

// 生成性能报告
void generate_performance_report(
    const std::vector<std::vector<PerformanceResult>> &all_results,
    const std::vector<u64> &test_sizes) {
  std::cout << "\n" << std::string(80, '=') << std::endl;
  std::cout << "                    GPU优化效果分析报告" << std::endl;
  std::cout << std::string(80, '=') << std::endl;

  // 按算法分类统计
  std::map<std::string, std::vector<double>> speedup_improvements;

  for (size_t i = 0; i < all_results.size(); i++) {
    u64 n = test_sizes[i];
    const auto &results = all_results[i];

    std::cout << "\n📊 n=" << n << " 性能分析:" << std::endl;

    // 按算法统计
    std::map<std::string, std::vector<PerformanceResult>> by_algorithm;
    for (const auto &result : results) {
      by_algorithm[result.algorithm_name].push_back(result);
    }

    for (const auto &[algo_name, algo_results] : by_algorithm) {
      if (algo_results.size() >= 2) { // 至少有原始和优化版本
        auto original_it =
            std::find_if(algo_results.begin(), algo_results.end(),
                         [](const PerformanceResult &r) {
                           return r.implementation_type == "GPU_Original";
                         });
        auto optimized_it =
            std::find_if(algo_results.begin(), algo_results.end(),
                         [](const PerformanceResult &r) {
                           return r.implementation_type == "GPU_Optimized";
                         });

        if (original_it != algo_results.end() &&
            optimized_it != algo_results.end()) {
          double improvement =
              original_it->execution_time_us / optimized_it->execution_time_us;
          speedup_improvements[algo_name].push_back(improvement);

          std::cout << "  " << algo_name << ": " << std::fixed
                    << std::setprecision(1) << original_it->execution_time_us
                    << "μs → " << optimized_it->execution_time_us << "μs"
                    << " (优化提升: " << std::setprecision(2) << improvement
                    << "x)" << std::endl;
        }
      }
    }
  }

  // 总体优化效果统计
  std::cout << "\n🚀 总体优化效果统计:" << std::endl;
  for (const auto &[algo_name, improvements] : speedup_improvements) {
    if (!improvements.empty()) {
      double avg_improvement =
          std::accumulate(improvements.begin(), improvements.end(), 0.0) /
          improvements.size();
      double max_improvement =
          *std::max_element(improvements.begin(), improvements.end());
      double min_improvement =
          *std::min_element(improvements.begin(), improvements.end());

      std::cout << "  " << algo_name << " 算法:" << std::endl;
      std::cout << "    平均优化提升: " << std::fixed << std::setprecision(2)
                << avg_improvement << "x" << std::endl;
      std::cout << "    最大优化提升: " << max_improvement << "x" << std::endl;
      std::cout << "    最小优化提升: " << min_improvement << "x" << std::endl;
    }
  }
}

int main() {
  std::cout << "=== GPU NTT优化效果对比测试 ===" << std::endl;

  // 检查CUDA设备
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0) {
    std::cout << "❌ 未发现CUDA设备!" << std::endl;
    return 1;
  }

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  std::cout << "🔧 使用GPU: " << prop.name << std::endl;
  std::cout << "🔧 显存: " << prop.totalGlobalMem / (1024 * 1024) << " MB"
            << std::endl;
  std::cout << "🔧 流式多处理器: " << prop.multiProcessorCount << std::endl;

  // 测试配置
  std::vector<u64> test_sizes = {1024, 4096, 16384, 65536};
  std::vector<std::string> algorithms = {"Montgomery", "Naive", "Barrett"};
  std::vector<std::vector<PerformanceResult>> all_results;

  // 测试参数
  u64 modulus = 7340033;
  u64 omega = 3;

  // 执行测试
  for (u64 n : test_sizes) {
    std::vector<PerformanceResult> size_results;

    for (const std::string &algo : algorithms) {
      auto algo_results =
          test_algorithm_comprehensive<u64>(n, modulus, omega, algo);
      size_results.insert(size_results.end(), algo_results.begin(),
                          algo_results.end());
    }

    all_results.push_back(size_results);
  }

  // 生成综合报告
  generate_performance_report(all_results, test_sizes);

  std::cout << "\n✅ 性能对比测试完成!" << std::endl;

  return 0;
}