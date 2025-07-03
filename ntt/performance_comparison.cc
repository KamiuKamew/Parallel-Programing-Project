#include "src/include/ntt_cuda.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <iomanip>
#include <functional>

using namespace std;
using namespace std::chrono;

// 高精度时间测量
double measureTime(std::function<void()> func)
{
    auto start = high_resolution_clock::now();
    func();
    auto end = high_resolution_clock::now();
    return duration_cast<nanoseconds>(end - start).count() / 1000.0; // 转换为微秒
}

// 生成测试数据
void generateTestData(u64 *a, u64 *b, u64 n, u64 p)
{
    for (u64 i = 0; i < n; i++)
    {
        a[i] = (i + 1) % p;
        b[i] = (i + 1) % p;
    }
}

// 验证结果一致性
bool verifyResults(u64 *result1, u64 *result2, u64 size)
{
    for (u64 i = 0; i < size; i++)
    {
        if (result1[i] != result2[i])
        {
            cout << "结果不一致：位置 " << i << " -> " << result1[i] << " vs " << result2[i] << endl;
            return false;
        }
    }
    return true;
}

int main()
{
    cout << "=== GPU NTT性能对比测试 ===" << endl;

    // 测试参数
    const u64 n = 131072;   // 大规模测试
    const u64 p = 7340033;  // 第一个测试模数
    const int num_runs = 5; // 多次运行取平均

    cout << "测试规模: n=" << n << ", p=" << p << endl;
    cout << "运行次数: " << num_runs << endl
         << endl;

    // 分配测试数据
    vector<u64> a(n), b(n);
    vector<u64> result_v1(2 * n - 1), result_v2_fixed(2 * n - 1), result_v2_enhanced(2 * n - 1);

    generateTestData(a.data(), b.data(), n, p);

    cout << "开始性能测试..." << endl;

    // 测试v1版本
    cout << "\n--- 测试v1版本（预计算旋转因子优化）---" << endl;
    double total_time_v1 = 0;
    for (int i = 0; i < num_runs; i++)
    {
        fill(result_v1.begin(), result_v1.end(), 0);
        double time = measureTime([&]()
                                  { poly_multiply_gpu_optimized_v1(a.data(), b.data(), result_v1.data(), n, p); });
        total_time_v1 += time;
        cout << "运行 " << (i + 1) << ": " << fixed << setprecision(2) << time << " μs" << endl;
    }
    double avg_time_v1 = total_time_v1 / num_runs;

    // 测试v2修复版本
    cout << "\n--- 测试v2修复版本（安全共享内存优化）---" << endl;
    double total_time_v2_fixed = 0;
    for (int i = 0; i < num_runs; i++)
    {
        fill(result_v2_fixed.begin(), result_v2_fixed.end(), 0);
        double time = measureTime([&]()
                                  { poly_multiply_gpu_optimized_v2_fixed(a.data(), b.data(), result_v2_fixed.data(), n, p); });
        total_time_v2_fixed += time;
        cout << "运行 " << (i + 1) << ": " << fixed << setprecision(2) << time << " μs" << endl;
    }
    double avg_time_v2_fixed = total_time_v2_fixed / num_runs;

    // 测试v2增强版本
    cout << "\n--- 测试v2增强版本（性能优化）---" << endl;
    double total_time_v2_enhanced = 0;
    for (int i = 0; i < num_runs; i++)
    {
        fill(result_v2_enhanced.begin(), result_v2_enhanced.end(), 0);
        double time = measureTime([&]()
                                  { poly_multiply_gpu_optimized_v2_enhanced(a.data(), b.data(), result_v2_enhanced.data(), n, p); });
        total_time_v2_enhanced += time;
        cout << "运行 " << (i + 1) << ": " << fixed << setprecision(2) << time << " μs" << endl;
    }
    double avg_time_v2_enhanced = total_time_v2_enhanced / num_runs;

    // 验证结果一致性
    cout << "\n=== 正确性验证 ===" << endl;
    bool v1_vs_v2_fixed = verifyResults(result_v1.data(), result_v2_fixed.data(), 2 * n - 1);
    bool v1_vs_v2_enhanced = verifyResults(result_v1.data(), result_v2_enhanced.data(), 2 * n - 1);
    bool v2_fixed_vs_v2_enhanced = verifyResults(result_v2_fixed.data(), result_v2_enhanced.data(), 2 * n - 1);

    cout << "v1 vs v2修复版: " << (v1_vs_v2_fixed ? "✅ 一致" : "❌ 不一致") << endl;
    cout << "v1 vs v2增强版: " << (v1_vs_v2_enhanced ? "✅ 一致" : "❌ 不一致") << endl;
    cout << "v2修复版 vs v2增强版: " << (v2_fixed_vs_v2_enhanced ? "✅ 一致" : "❌ 不一致") << endl;

    // 性能总结
    cout << "\n=== 性能总结 ===" << endl;
    cout << fixed << setprecision(2);
    cout << "v1版本平均时间:        " << avg_time_v1 << " μs" << endl;
    cout << "v2修复版平均时间:      " << avg_time_v2_fixed << " μs" << endl;
    cout << "v2增强版平均时间:      " << avg_time_v2_enhanced << " μs" << endl;

    cout << "\n=== 性能提升分析 ===" << endl;
    double speedup_v2_fixed = avg_time_v1 / avg_time_v2_fixed;
    double speedup_v2_enhanced = avg_time_v1 / avg_time_v2_enhanced;
    double improvement_enhanced = avg_time_v2_fixed / avg_time_v2_enhanced;

    cout << "v2修复版相对v1提升:    " << speedup_v2_fixed << "x" << endl;
    cout << "v2增强版相对v1提升:    " << speedup_v2_enhanced << "x" << endl;
    cout << "v2增强版相对v2修复版:  " << improvement_enhanced << "x" << endl;

    // 推荐建议
    cout << "\n=== 推荐建议 ===" << endl;
    if (avg_time_v2_enhanced < avg_time_v2_fixed && v1_vs_v2_enhanced)
    {
        cout << "🏆 推荐使用v2增强版：最佳性能且正确性验证通过" << endl;
    }
    else if (v1_vs_v2_fixed)
    {
        cout << "🛡️ 推荐使用v2修复版：稳定可靠的性能提升" << endl;
    }
    else
    {
        cout << "⚡ 推荐使用v1版本：最可靠的基准实现" << endl;
    }

    return 0;
}