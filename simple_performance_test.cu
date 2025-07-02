#include "ntt/src/include/ntt.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std::chrono;

// GPU优化版本函数声明（包含在最后）
template <typename T>
void poly_multiply_ntt_gpu_mont_optimized(T *a, T *b, T *ab, T n, T p,
                                          T omega = 3);

struct BenchmarkResult {
  uint32_t n;
  double cpu_time_us;
  double gpu_opt_time_us;
  double speedup;
  bool correctness_pass;
};

template <typename T> void generate_test_data(T *a, T *b, T n) {
  for (T i = 0; i < n; i++) {
    a[i] = (i * 17 + 7) % 1000 + 1;
    b[i] = (i * 23 + 11) % 1000 + 1;
  }
}

template <typename T>
bool verify_correctness(T *cpu_result, T *gpu_result, T result_size) {
  for (T i = 0; i < result_size; i++) {
    if (cpu_result[i] != gpu_result[i]) {
      return false;
    }
  }
  return true;
}

void run_benchmark() {
  const uint32_t p = 998244353;
  const uint32_t omega = 3;

  std::vector<uint32_t> test_sizes = {1024, 4096, 16384, 65536, 131072, 262144};
  std::vector<BenchmarkResult> results;

  std::cout << "=== GPU NTT 优化性能验证 ===" << std::endl;
  std::cout << "对比：CPU基线版本 vs GPU优化版本" << std::endl;
  std::cout << std::string(70, '=') << std::endl;
  std::cout << std::setw(10) << "Size" << std::setw(15) << "CPU (μs)"
            << std::setw(15) << "GPU优化 (μs)" << std::setw(12) << "加速比"
            << std::setw(10) << "正确性" << std::endl;
  std::cout << std::string(70, '-') << std::endl;

  for (uint32_t n : test_sizes) {
    BenchmarkResult result;
    result.n = n;

    // 分配内存
    uint32_t *a = new uint32_t[n];
    uint32_t *b = new uint32_t[n];
    uint32_t *cpu_result = new uint32_t[2 * n - 1]();
    uint32_t *gpu_result = new uint32_t[2 * n - 1]();

    // 生成测试数据
    generate_test_data(a, b, n);

    // CPU版本测试 - 运行多次取平均
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

    // GPU优化版本预热
    poly_multiply_ntt_gpu_mont_optimized(a, b, gpu_result, n, p, omega);

    // GPU优化版本测试
    const int gpu_runs = (n <= 16384) ? 3 : 1;
    double gpu_total_time = 0.0;

    for (int run = 0; run < gpu_runs; run++) {
      uint32_t *a_copy = new uint32_t[n];
      uint32_t *b_copy = new uint32_t[n];
      std::copy(a, a + n, a_copy);
      std::copy(b, b + n, b_copy);

      auto gpu_start = high_resolution_clock::now();
      poly_multiply_ntt_gpu_mont_optimized(a_copy, b_copy, gpu_result, n, p,
                                           omega);
      auto gpu_end = high_resolution_clock::now();

      gpu_total_time +=
          duration_cast<microseconds>(gpu_end - gpu_start).count();

      delete[] a_copy;
      delete[] b_copy;
    }
    result.gpu_opt_time_us = gpu_total_time / gpu_runs;

    // 验证正确性
    result.correctness_pass =
        verify_correctness(cpu_result, gpu_result, 2 * n - 1);
    result.speedup = result.cpu_time_us / result.gpu_opt_time_us;

    results.push_back(result);

    // 输出结果
    std::cout << std::setw(10) << result.n << std::setw(15) << std::fixed
              << std::setprecision(1) << result.cpu_time_us << std::setw(15)
              << std::fixed << std::setprecision(1) << result.gpu_opt_time_us
              << std::setw(12) << std::fixed << std::setprecision(2)
              << result.speedup << std::setw(10)
              << (result.correctness_pass ? "✓" : "✗") << std::endl;

    delete[] a;
    delete[] b;
    delete[] cpu_result;
    delete[] gpu_result;
  }

  // 性能统计
  std::cout << std::string(70, '=') << std::endl;
  std::cout << "=== 优化效果分析 ===" << std::endl;

  double total_speedup = 0.0;
  int valid_count = 0;
  double max_speedup = 0.0;
  uint32_t max_speedup_size = 0;

  // 分析大规模问题性能
  double large_scale_speedup = 0.0;
  int large_count = 0;

  for (const auto &res : results) {
    if (res.correctness_pass) {
      total_speedup += res.speedup;
      valid_count++;

      if (res.speedup > max_speedup) {
        max_speedup = res.speedup;
        max_speedup_size = res.n;
      }

      if (res.n >= 65536) {
        large_scale_speedup += res.speedup;
        large_count++;
      }
    }
  }

  if (valid_count > 0) {
    std::cout << "📊 整体平均加速比: " << std::fixed << std::setprecision(2)
              << (total_speedup / valid_count) << "x" << std::endl;
    std::cout << "🚀 最大加速比: " << std::fixed << std::setprecision(2)
              << max_speedup << "x (n=" << max_speedup_size << ")" << std::endl;

    if (large_count > 0) {
      std::cout << "📈 大规模问题平均加速比 (n≥65536): " << std::fixed
                << std::setprecision(2) << (large_scale_speedup / large_count)
                << "x" << std::endl;
    }

    // 评估优化效果
    double avg_speedup = total_speedup / valid_count;
    std::cout << std::endl << "🎯 优化效果评估: ";
    if (avg_speedup >= 3.0) {
      std::cout << "显著优化 - GPU充分发挥并行优势!" << std::endl;
    } else if (avg_speedup >= 2.0) {
      std::cout << "良好优化 - GPU实现了预期加速效果" << std::endl;
    } else if (avg_speedup >= 1.5) {
      std::cout << "一般优化 - 仍有进一步提升空间" << std::endl;
    } else {
      std::cout << "优化不足 - 需要深入分析瓶颈" << std::endl;
    }
  }

  std::cout << std::endl << "CSV数据 (用于绘图分析):" << std::endl;
  std::cout << "Size,CPU_Time_us,GPU_Opt_Time_us,Speedup,Correct" << std::endl;
  for (const auto &res : results) {
    std::cout << res.n << "," << res.cpu_time_us << "," << res.gpu_opt_time_us
              << "," << res.speedup << ","
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
  std::cout << std::endl;

  try {
    run_benchmark();
  } catch (const std::exception &e) {
    std::cerr << "❌ 测试错误: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}

// 最后包含GPU优化实现
#include "ntt/src/include/CUDA/ntt_gpu_optimized.cu"