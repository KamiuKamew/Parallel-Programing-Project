#include "src/include/PerformanceAnalysis/comprehensive_analysis.h"
#include "src/include/general/utils.h"
#include <iostream>
#include <vector>

using u64 = uint64_t;

int main()
{
    std::cout << "=== 性能分析测试程序 ===" << std::endl;

    // 测试参数
    const u64 n = 1024;      // 较小的问题规模以便快速测试
    const u64 p = 998244353; // NTT友好的素数
    const u64 omega = 3;     // 原根
    const std::vector<int> thread_counts = {1, 2, 4};

    std::cout << "\n1. 测试线程创建开销分析..." << std::endl;
    for (int threads = 2; threads <= 4; threads += 2)
    {
        ThreadOverheadResult result = measure_thread_overhead<u64>(n, p, omega, threads, 3);
        std::cout << "  " << threads << "线程: 线程池=" << result.thread_pool_time_ms
                  << "ms, 重复创建=" << result.thread_creation_time_ms
                  << "ms, 开销比例=" << result.overhead_ratio << "x" << std::endl;
    }

    std::cout << "\n2. 测试同步开销分析..." << std::endl;
    auto sync_results = analyze_sync_overhead_scaling<u64>(n, p, omega, thread_counts, 2);
    for (const auto &result : sync_results)
    {
        std::cout << "  " << result.num_threads << "线程: 总时间=" << result.total_time_ms
                  << "ms, 同步时间=" << result.sync_time_ms
                  << "ms, 同步占比=" << (result.sync_overhead_ratio * 100) << "%" << std::endl;
    }

    std::cout << "\n3. 测试内存带宽分析..." << std::endl;
    auto memory_analysis = analyze_memory_contention(thread_counts, 64, 100); // 很小的参数以加快测试
    std::cout << "  顺序访问结果:" << std::endl;
    for (const auto &result : memory_analysis.sequential_results)
    {
        std::cout << "    " << result.num_threads << "线程: " << result.bandwidth_gbps
                  << " GB/s, 效率=" << (result.efficiency_ratio * 100) << "%" << std::endl;
    }

    std::cout << "\n4. 测试False Sharing分析..." << std::endl;
    for (int threads = 2; threads <= 4; threads += 2)
    {
        CacheFriendlinessResult result = test_false_sharing_impact(threads, 100000); // 较小迭代数
        std::cout << "  " << threads << "线程: 对齐=" << result.aligned_time_ms
                  << "ms, 未对齐=" << result.unaligned_time_ms
                  << "ms, false sharing开销=" << result.false_sharing_overhead << "x" << std::endl;
    }

    std::cout << "\n5. 测试调度策略对比..." << std::endl;
    auto scheduling_result = compare_scheduling_overhead<u64>(n, p, omega, 4);
    std::cout << "  静态调度: " << scheduling_result.static_result.total_time_ms << "ms" << std::endl;
    std::cout << "  动态调度: " << scheduling_result.dynamic_result.total_time_ms << "ms (开销比例: "
              << scheduling_result.dynamic_overhead_ratio << "x)" << std::endl;
    std::cout << "  guided调度: " << scheduling_result.guided_result.total_time_ms << "ms (开销比例: "
              << scheduling_result.guided_overhead_ratio << "x)" << std::endl;

    // 测试CSV导出功能
    std::cout << "\n6. 生成CSV报告..." << std::endl;
    ComprehensivePerformanceReport report;

    // 填充一些测试数据
    for (int threads = 2; threads <= 4; threads += 2)
    {
        ThreadOverheadResult thread_result = measure_thread_overhead<u64>(n, p, omega, threads, 2);
        report.thread_overhead_results.push_back(thread_result);
    }

    report.sync_overhead_results = analyze_sync_overhead_scaling<u64>(n, p, omega, thread_counts, 2);
    report.memory_analysis = analyze_memory_contention(thread_counts, 32, 50);

    export_detailed_csv_report(report, "performance_analysis_test.csv");

    std::cout << "\n=== 所有测试完成 ===" << std::endl;

    return 0;
}