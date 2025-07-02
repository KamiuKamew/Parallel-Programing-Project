#include "src/include/CUDA/ntt.h"
#include "src/include/ntt.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

int main() {
  cout << "=========================================" << endl;
  cout << "三种模乘算法全面测试" << endl;
  cout << "=========================================" << endl;

  // 测试参数
  u64 n = 4;
  u64 p = 998244353;
  u64 omega = 3;

  // 测试数据
  u64 a[4] = {1, 2, 3, 4};
  u64 b[4] = {1, 1, 1, 1};
  u64 ab_cpu[7], ab_mont[7], ab_naive[7], ab_barrett[7];

  cout << "测试数据: a=[1,2,3,4], b=[1,1,1,1]" << endl;
  cout << "预期结果: [1,3,6,10,9,7,4]" << endl;

  // 1. CPU基线测试
  cout << "\n=== CPU基线测试 ===" << endl;
  auto start_cpu = chrono::high_resolution_clock::now();
  poly_multiply_ntt(a, b, ab_cpu, n, p, omega);
  auto end_cpu = chrono::high_resolution_clock::now();
  auto cpu_time =
      chrono::duration_cast<chrono::microseconds>(end_cpu - start_cpu).count();

  cout << "CPU结果: [";
  for (int i = 0; i < 7; i++)
    cout << ab_cpu[i] << " ";
  cout << "]" << endl;
  cout << "CPU用时: " << cpu_time << "μs" << endl;

  // 2. GPU Montgomery测试
  cout << "\n=== GPU Montgomery模乘测试 ===" << endl;
  auto start_mont = chrono::high_resolution_clock::now();
  poly_multiply_ntt_gpu_mont(a, b, ab_mont, n, p, omega);
  auto end_mont = chrono::high_resolution_clock::now();
  auto mont_time =
      chrono::duration_cast<chrono::microseconds>(end_mont - start_mont)
          .count();

  cout << "GPU Montgomery: [";
  for (int i = 0; i < 7; i++)
    cout << ab_mont[i] << " ";
  cout << "]" << endl;
  cout << "Montgomery用时: " << mont_time << "μs" << endl;

  bool mont_correct = true;
  for (int i = 0; i < 7; i++) {
    if (ab_mont[i] != ab_cpu[i]) {
      mont_correct = false;
      break;
    }
  }
  cout << "Montgomery正确性: " << (mont_correct ? "✅" : "❌") << endl;

  // 3. GPU Naive测试
  cout << "\n=== GPU Naive模乘测试 ===" << endl;
  auto start_naive = chrono::high_resolution_clock::now();
  poly_multiply_ntt_gpu_naive(a, b, ab_naive, n, p, omega);
  auto end_naive = chrono::high_resolution_clock::now();
  auto naive_time =
      chrono::duration_cast<chrono::microseconds>(end_naive - start_naive)
          .count();

  cout << "GPU Naive: [";
  for (int i = 0; i < 7; i++)
    cout << ab_naive[i] << " ";
  cout << "]" << endl;
  cout << "Naive用时: " << naive_time << "μs" << endl;

  bool naive_correct = true;
  for (int i = 0; i < 7; i++) {
    if (ab_naive[i] != ab_cpu[i]) {
      naive_correct = false;
      break;
    }
  }
  cout << "Naive正确性: " << (naive_correct ? "✅" : "❌") << endl;

  // 4. GPU Barrett测试
  cout << "\n=== GPU Barrett模乘测试 ===" << endl;
  auto start_barrett = chrono::high_resolution_clock::now();
  poly_multiply_ntt_gpu_barrett(a, b, ab_barrett, n, p, omega);
  auto end_barrett = chrono::high_resolution_clock::now();
  auto barrett_time =
      chrono::duration_cast<chrono::microseconds>(end_barrett - start_barrett)
          .count();

  cout << "GPU Barrett: [";
  for (int i = 0; i < 7; i++)
    cout << ab_barrett[i] << " ";
  cout << "]" << endl;
  cout << "Barrett用时: " << barrett_time << "μs" << endl;

  bool barrett_correct = true;
  for (int i = 0; i < 7; i++) {
    if (ab_barrett[i] != ab_cpu[i]) {
      barrett_correct = false;
      break;
    }
  }
  cout << "Barrett正确性: " << (barrett_correct ? "✅" : "❌") << endl;

  // 性能对比分析
  cout << "\n=========================================" << endl;
  cout << "三种模乘算法性能对比分析" << endl;
  cout << "=========================================" << endl;

  cout << left << setw(15) << "算法" << setw(12) << "时间(μs)" << setw(12)
       << "加速比" << setw(10) << "正确性" << endl;
  cout << string(50, '-') << endl;

  cout << left << setw(15) << "CPU Montgomery" << setw(12) << cpu_time
       << setw(12) << "1.000x" << setw(10) << "✅" << endl;

  double mont_speedup = (double)cpu_time / mont_time;
  cout << left << setw(15) << "GPU Montgomery" << setw(12) << mont_time
       << setw(12) << fixed << setprecision(3) << mont_speedup << "x"
       << setw(10) << (mont_correct ? "✅" : "❌") << endl;

  double naive_speedup = (double)cpu_time / naive_time;
  cout << left << setw(15) << "GPU Naive" << setw(12) << naive_time << setw(12)
       << fixed << setprecision(3) << naive_speedup << "x" << setw(10)
       << (naive_correct ? "✅" : "❌") << endl;

  double barrett_speedup = (double)cpu_time / barrett_time;
  cout << left << setw(15) << "GPU Barrett" << setw(12) << barrett_time
       << setw(12) << fixed << setprecision(3) << barrett_speedup << "x"
       << setw(10) << (barrett_correct ? "✅" : "❌") << endl;

  // 关键发现
  cout << "\n=== 关键发现和分析 ===" << endl;

  cout << "1. 功能正确性验证:" << endl;
  cout << "   - Montgomery模乘: "
       << (mont_correct ? "✅ 完全正确" : "❌ 存在错误") << endl;
  cout << "   - Naive模乘:      "
       << (naive_correct ? "✅ 完全正确" : "❌ 存在错误") << endl;
  cout << "   - Barrett模乘:    "
       << (barrett_correct ? "✅ 完全正确" : "❌ 存在错误") << endl;

  cout << "\n2. 性能特征分析:" << endl;
  cout << "   - CPU基线性能: " << cpu_time << "μs" << endl;
  cout << "   - 最快GPU算法: ";

  vector<pair<string, long long>> gpu_times = {{"Montgomery", mont_time},
                                               {"Naive", naive_time},
                                               {"Barrett", barrett_time}};

  sort(gpu_times.begin(), gpu_times.end(),
       [](const auto &a, const auto &b) { return a.second < b.second; });

  cout << gpu_times[0].first << " (" << gpu_times[0].second << "μs)" << endl;

  cout << "\n3. 规约算法验证:" << endl;
  if (barrett_correct && barrett_time < mont_time) {
    cout << "   ✅ Barrett规约算法比基础Montgomery模乘更快" << endl;
    cout << "   ✅ 性能提升: " << fixed << setprecision(2)
         << (double)mont_time / barrett_time << "x" << endl;
  } else if (barrett_correct) {
    cout << "   ⚠️  Barrett规约算法功能正确，但性能优势不明显" << endl;
  } else {
    cout << "   ❌ Barrett规约算法存在正确性问题" << endl;
  }

  cout << "\n4. 实验要求验证:" << endl;
  cout << "   - ✅ 对比了基础模乘和两种优化模乘在GPU上的加速比" << endl;
  cout << "   - " << (barrett_time < mont_time ? "✅" : "❌")
       << " 验证了规约算法对模乘优化的加速效果" << endl;

  cout << "\n=========================================" << endl;
  cout << "三种模乘算法测试完成" << endl;
  cout << "所有算法正确性: "
       << ((mont_correct && naive_correct && barrett_correct) ? "✅ 全部正确"
                                                              : "❌ 存在错误")
       << endl;
  cout << "=========================================" << endl;

  return (mont_correct && naive_correct && barrett_correct) ? 0 : 1;
}