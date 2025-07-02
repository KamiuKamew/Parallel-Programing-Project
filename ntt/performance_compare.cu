#include "src/include/CUDA/ntt.h"
#include "src/include/ntt.h"
#include "src/include/transform.h"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

// 从文件读取测试数据
bool load_test_data(const std::string &filename, std::vector<u64> &a,
                    std::vector<u64> &b, u64 &n, u64 &p) {
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    printf("无法打开文件: %s\n", filename.c_str());
    return false;
  }

  fin >> n >> p;
  printf("数据规模: n=%lu, p=%lu\n", n, p);

  a.resize(n);
  b.resize(n);

  for (u64 i = 0; i < n; i++) {
    fin >> a[i];
  }

  for (u64 i = 0; i < n; i++) {
    fin >> b[i];
  }

  fin.close();
  return true;
}

// 验证结果是否一致
bool verify_results(const std::vector<u64> &result1,
                    const std::vector<u64> &result2, u64 result_size) {
  for (u64 i = 0; i < result_size; i++) {
    if (result1[i] != result2[i]) {
      printf("结果不匹配: 位置%lu, 值1=%lu, 值2=%lu\n", i, result1[i],
             result2[i]);
      return false;
    }
  }
  return true;
}

int main() {
  printf("=== GPU vs CPU NTT 性能对比测试 ===\n\n");

  // 检查CUDA设备
  int device_count;
  cudaError_t err = cudaGetDeviceCount(&device_count);
  if (err != cudaSuccess || device_count == 0) {
    printf("错误: 未检测到CUDA设备\n");
    return 1;
  }
  printf("检测到 %d 个CUDA设备\n\n", device_count);

  // 加载测试数据 - 使用1.in (n=131072)
  std::vector<u64> a, b;
  u64 n, p;
  std::string test_file =
      "/home/hexay/projects/Lab Proj/Parallel-Programing-Project/.nttdata/1.in";

  if (!load_test_data(test_file, a, b, n, p)) {
    return 1;
  }

  // 原根设置 (根据main.cc的注释，所有模数的原根都是3)
  u64 omega = 3;

  // 准备结果数组
  u64 result_size = 2 * n - 1;
  std::vector<u64> result_cpu(result_size, 0);
  std::vector<u64> result_gpu(result_size, 0);

  // 为了公平对比，使用相同的输入数据副本
  std::vector<u64> a_cpu = a;
  std::vector<u64> b_cpu = b;
  std::vector<u64> a_gpu = a;
  std::vector<u64> b_gpu = b;

  printf("开始性能测试...\n\n");

  // CPU版本性能测试
  printf("--- CPU版本测试 ---\n");
  auto start = std::chrono::high_resolution_clock::now();

  poly_multiply_ntt(a_cpu.data(), b_cpu.data(), result_cpu.data(), n, p, omega);

  auto end = std::chrono::high_resolution_clock::now();
  auto cpu_time =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count();

  printf("CPU时间: %ld μs\n\n", cpu_time);

  // GPU版本性能测试
  printf("--- GPU版本测试 ---\n");
  start = std::chrono::high_resolution_clock::now();

  poly_multiply_ntt_gpu(a_gpu.data(), b_gpu.data(), result_gpu.data(), n, p,
                        omega);

  end = std::chrono::high_resolution_clock::now();
  auto gpu_time =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count();

  printf("GPU时间: %ld μs\n\n", gpu_time);

  // 验证结果正确性
  printf("--- 结果验证 ---\n");
  bool results_match = verify_results(result_cpu, result_gpu, result_size);

  if (results_match) {
    printf("✓ GPU和CPU结果一致\n");
  } else {
    printf("✗ GPU和CPU结果不一致\n");
    return 1;
  }

  // 性能总结
  printf("\n=== 性能总结 ===\n");
  printf("数据规模: n = %lu\n", n);
  printf("CPU时间:  %ld μs\n", cpu_time);
  printf("GPU时间:  %ld μs\n", gpu_time);

  if (gpu_time > 0) {
    double speedup = (double)cpu_time / gpu_time;
    printf("加速比:   %.2fx", speedup);
    if (speedup > 1.0) {
      printf(" (GPU更快)\n");
    } else {
      printf(" (CPU更快)\n");
    }
  }

  printf("GPU开销分析:\n");
  printf("- 数据传输开销包含在总时间内\n");
  printf("- 大规模数据下GPU并行优势更明显\n");

  return 0;
}