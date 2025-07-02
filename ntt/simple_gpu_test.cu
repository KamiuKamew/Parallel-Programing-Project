#include "src/include/CUDA/ntt.h"
#include "src/include/ntt.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>

// 简单的计时器
class SimpleTimer {
private:
  std::chrono::high_resolution_clock::time_point start_time;

public:
  void start() { start_time = std::chrono::high_resolution_clock::now(); }

  double stop() {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
        end_time - start_time);
    return duration.count();
  }
};

// 测试单个算法
void test_algorithm(const std::string &name,
                    void (*gpu_func)(u64 *, u64 *, u64 *, u64, u64, u64), u64 n,
                    u64 p, u64 omega) {
  std::cout << "\n=== 测试 " << name << " 算法 (n=" << n
            << ") ===" << std::endl;

  // 生成测试数据
  u64 *a = new u64[n];
  u64 *b = new u64[n];
  u64 *cpu_result = new u64[2 * n - 1];
  u64 *gpu_result = new u64[2 * n - 1];

  for (u64 i = 0; i < n; i++) {
    a[i] = (i * 17 + 7) % 1000 + 1;
    b[i] = (i * 23 + 11) % 1000 + 1;
  }

  // 初始化结果数组
  for (u64 i = 0; i < 2 * n - 1; i++) {
    cpu_result[i] = 0;
    gpu_result[i] = 0;
  }

  SimpleTimer timer;

  try {
    // CPU版本测试（只对Montgomery有效）
    double cpu_time = -1;
    if (name == "Montgomery") {
      u64 *a_copy = new u64[n];
      u64 *b_copy = new u64[n];
      std::copy(a, a + n, a_copy);
      std::copy(b, b + n, b_copy);

      timer.start();
      poly_multiply_ntt(a_copy, b_copy, cpu_result, n, p, omega);
      cpu_time = timer.stop();

      delete[] a_copy;
      delete[] b_copy;

      std::cout << "CPU时间: " << cpu_time << " μs" << std::endl;
    }

    // GPU版本测试
    u64 *a_copy = new u64[n];
    u64 *b_copy = new u64[n];
    std::copy(a, a + n, a_copy);
    std::copy(b, b + n, b_copy);

    timer.start();
    gpu_func(a_copy, b_copy, gpu_result, n, p, omega);
    double gpu_time = timer.stop();

    delete[] a_copy;
    delete[] b_copy;

    std::cout << "GPU时间: " << gpu_time << " μs" << std::endl;

    // 计算加速比
    if (cpu_time > 0) {
      double speedup = cpu_time / gpu_time;
      std::cout << "加速比: " << speedup << "x" << std::endl;
    }

    // 验证正确性（只对Montgomery有效）
    if (name == "Montgomery") {
      bool correct = true;
      for (u64 i = 0; i < 2 * n - 1; i++) {
        if (cpu_result[i] != gpu_result[i]) {
          correct = false;
          std::cout << "结果不匹配 at index " << i << ": CPU=" << cpu_result[i]
                    << " GPU=" << gpu_result[i] << std::endl;
          break;
        }
      }

      if (correct) {
        std::cout << "✅ 正确性验证通过" << std::endl;
      } else {
        std::cout << "❌ 正确性验证失败" << std::endl;
      }
    } else {
      std::cout << "✅ GPU执行成功" << std::endl;
    }

  } catch (const std::exception &e) {
    std::cout << "❌ 异常: " << e.what() << std::endl;
  } catch (...) {
    std::cout << "❌ 未知异常" << std::endl;
  }

  delete[] a;
  delete[] b;
  delete[] cpu_result;
  delete[] gpu_result;
}

int main() {
  std::cout << "=== 简化GPU NTT测试 ===" << std::endl;

  // 检查CUDA设备
  int deviceCount;
  if (cudaGetDeviceCount(&deviceCount) != cudaSuccess || deviceCount == 0) {
    std::cout << "❌ 未发现CUDA设备!" << std::endl;
    return 1;
  }

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  std::cout << "🔧 GPU: " << prop.name << std::endl;
  std::cout << "🔧 显存: " << prop.totalGlobalMem / (1024 * 1024) << " MB"
            << std::endl;

  // 测试参数
  u64 modulus = 7340033;
  u64 omega = 3;

  // 小规模测试
  std::vector<u64> test_sizes = {1024, 4096};

  for (u64 n : test_sizes) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "问题规模: n=" << n << std::endl;

    // 测试Montgomery算法
    test_algorithm("Montgomery", poly_multiply_ntt_gpu_mont, n, modulus, omega);

    // 测试Naive算法
    test_algorithm("Naive", poly_multiply_ntt_gpu_naive, n, modulus, omega);

    // 测试Barrett算法
    test_algorithm("Barrett", poly_multiply_ntt_gpu_barrett, n, modulus, omega);
  }

  std::cout << "\n✅ 测试完成!" << std::endl;
  return 0;
}