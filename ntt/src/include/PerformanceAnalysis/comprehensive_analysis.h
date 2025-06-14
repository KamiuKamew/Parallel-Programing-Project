#pragma once

#include "thread_overhead.h"
#include "sync_overhead.h"
#include "memory_bandwidth.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

/**
 * @brief 综合性能分析结果
 */
struct ComprehensivePerformanceReport
{
    std::vector<ThreadOverheadResult> thread_overhead_results;
    std::vector<SyncOverheadResult> sync_overhead_results;
    MemoryContentionAnalysis memory_analysis;
    std::vector<CacheFriendlinessResult> cache_results;

    // 关键发现
    int optimal_thread_count_for_computation;
    int optimal_thread_count_for_memory;
    double max_thread_creation_overhead;
    double max_sync_overhead_ratio;
    double max_false_sharing_overhead;

    // 性能退化分析
    struct PerformanceDegradationAnalysis
    {
        bool thread_creation_is_bottleneck;
        bool sync_overhead_is_bottleneck;
        bool memory_bandwidth_is_bottleneck;
        bool false_sharing_is_bottleneck;
        std::string primary_bottleneck;
        std::vector<std::string> recommendations;
    } degradation_analysis;
};

/**
 * @brief 执行全面性能分析
 */
template <typename T>
ComprehensivePerformanceReport run_comprehensive_analysis(T n, T p, T omega,
                                                          const std::vector<int> &thread_counts,
                                                          bool enable_detailed_output = true)
{
    ComprehensivePerformanceReport report;

    if (enable_detailed_output)
    {
        std::cout << "=== 开始综合性能分析 ===" << std::endl;
        std::cout << "问题规模: n=" << n << ", 模数: p=" << p << std::endl;
        std::cout << "测试线程数: ";
        for (int t : thread_counts)
            std::cout << t << " ";
        std::cout << std::endl
                  << std::endl;
    }

    // 1. 线程创建开销分析
    if (enable_detailed_output)
    {
        std::cout << "1. 分析线程创建开销..." << std::endl;
    }

    for (int threads : thread_counts)
    {
        if (threads == 1)
            continue; // 单线程不需要测试线程创建开销
        ThreadOverheadResult result = measure_thread_overhead<T>(n, p, omega, threads, 5);
        report.thread_overhead_results.push_back(result);

        if (enable_detailed_output)
        {
            std::cout << "  " << threads << "线程: 线程池=" << std::fixed << std::setprecision(2)
                      << result.thread_pool_time_ms << "ms, 重复创建=" << result.thread_creation_time_ms
                      << "ms, 开销比例=" << result.overhead_ratio << "x" << std::endl;
        }
    }

    // 2. 同步开销分析
    if (enable_detailed_output)
    {
        std::cout << "\n2. 分析同步开销..." << std::endl;
    }

    report.sync_overhead_results = analyze_sync_overhead_scaling<T>(n, p, omega, thread_counts, 3);

    for (const auto &result : report.sync_overhead_results)
    {
        if (enable_detailed_output)
        {
            std::cout << "  " << result.num_threads << "线程: 总时间=" << std::fixed << std::setprecision(2)
                      << result.total_time_ms << "ms, 同步时间=" << result.sync_time_ms
                      << "ms, 同步占比=" << (result.sync_overhead_ratio * 100) << "%" << std::endl;
        }
    }

    // 3. 内存带宽竞争分析
    if (enable_detailed_output)
    {
        std::cout << "\n3. 分析内存带宽竞争..." << std::endl;
    }

    report.memory_analysis = analyze_memory_contention(thread_counts, 128, 500); // 较小的内存和迭代数以加快测试

    if (enable_detailed_output)
    {
        std::cout << "  顺序访问结果:" << std::endl;
        for (const auto &result : report.memory_analysis.sequential_results)
        {
            std::cout << "    " << result.num_threads << "线程: " << std::fixed << std::setprecision(2)
                      << result.bandwidth_gbps << " GB/s, 效率=" << (result.efficiency_ratio * 100) << "%" << std::endl;
        }
    }

    // 4. Cache友好性分析
    if (enable_detailed_output)
    {
        std::cout << "\n4. 分析False Sharing影响..." << std::endl;
    }

    for (int threads : thread_counts)
    {
        if (threads == 1)
            continue;
        CacheFriendlinessResult result = test_false_sharing_impact(threads, 1000000); // 较小迭代数
        report.cache_results.push_back(result);

        if (enable_detailed_output)
        {
            std::cout << "  " << threads << "线程: 对齐=" << std::fixed << std::setprecision(2)
                      << result.aligned_time_ms << "ms, 未对齐=" << result.unaligned_time_ms
                      << "ms, false sharing开销=" << result.false_sharing_overhead << "x" << std::endl;
        }
    }

    // 5. 综合分析和建议生成
    // analyze_performance_bottlenecks(report);  // 暂时注释掉，稍后在函数定义后调用

    // if (enable_detailed_output)
    // {
    //     print_comprehensive_summary(report);  // 暂时注释掉，稍后在函数定义后调用
    // }

    return report;
}

/**
 * @brief 分析性能瓶颈
 */
void analyze_performance_bottlenecks(ComprehensivePerformanceReport &report)
{
    auto &analysis = report.degradation_analysis;

    // 找出最大的各类开销
    report.max_thread_creation_overhead = 0.0;
    for (const auto &result : report.thread_overhead_results)
    {
        report.max_thread_creation_overhead = std::max(report.max_thread_creation_overhead, result.overhead_ratio);
    }

    report.max_sync_overhead_ratio = 0.0;
    for (const auto &result : report.sync_overhead_results)
    {
        report.max_sync_overhead_ratio = std::max(report.max_sync_overhead_ratio, result.sync_overhead_ratio);
    }

    report.max_false_sharing_overhead = 0.0;
    for (const auto &result : report.cache_results)
    {
        report.max_false_sharing_overhead = std::max(report.max_false_sharing_overhead, result.false_sharing_overhead);
    }

    // 判断主要瓶颈
    const double THREAD_CREATION_THRESHOLD = 1.5; // 1.5x以上认为是显著开销
    const double SYNC_OVERHEAD_THRESHOLD = 0.2;   // 20%以上认为是显著开销
    const double FALSE_SHARING_THRESHOLD = 1.3;   // 1.3x以上认为是显著开销

    analysis.thread_creation_is_bottleneck = (report.max_thread_creation_overhead > THREAD_CREATION_THRESHOLD);
    analysis.sync_overhead_is_bottleneck = (report.max_sync_overhead_ratio > SYNC_OVERHEAD_THRESHOLD);
    analysis.false_sharing_is_bottleneck = (report.max_false_sharing_overhead > FALSE_SHARING_THRESHOLD);

    // 内存带宽瓶颈判断 - 如果多线程效率低于80%
    analysis.memory_bandwidth_is_bottleneck = false;
    for (const auto &result : report.memory_analysis.sequential_results)
    {
        if (result.num_threads > 1 && result.efficiency_ratio < 0.8)
        {
            analysis.memory_bandwidth_is_bottleneck = true;
            break;
        }
    }

    // 确定主要瓶颈
    if (analysis.thread_creation_is_bottleneck)
    {
        analysis.primary_bottleneck = "线程创建开销";
    }
    else if (analysis.sync_overhead_is_bottleneck)
    {
        analysis.primary_bottleneck = "线程同步开销";
    }
    else if (analysis.memory_bandwidth_is_bottleneck)
    {
        analysis.primary_bottleneck = "内存带宽竞争";
    }
    else if (analysis.false_sharing_is_bottleneck)
    {
        analysis.primary_bottleneck = "False Sharing";
    }
    else
    {
        analysis.primary_bottleneck = "其他因素";
    }

    // 生成优化建议
    if (analysis.thread_creation_is_bottleneck)
    {
        analysis.recommendations.push_back("使用线程池复用线程，避免频繁创建/销毁");
        analysis.recommendations.push_back("减少并行区域的数量，增加每个区域的工作量");
    }
    if (analysis.sync_overhead_is_bottleneck)
    {
        analysis.recommendations.push_back("减少同步点的数量");
        analysis.recommendations.push_back("使用静态调度替代动态调度");
        analysis.recommendations.push_back("增加每个线程的工作粒度");
    }
    if (analysis.memory_bandwidth_is_bottleneck)
    {
        analysis.recommendations.push_back("减少线程数以避免内存带宽饱和");
        analysis.recommendations.push_back("使用NUMA感知的内存分配");
        analysis.recommendations.push_back("优化数据访问模式以提高缓存命中率");
    }
    if (analysis.false_sharing_is_bottleneck)
    {
        analysis.recommendations.push_back("使用缓存行对齐(alignas(64))避免false sharing");
        analysis.recommendations.push_back("重新设计数据结构以减少共享数据");
    }
}

/**
 * @brief 打印综合分析摘要
 */
void print_comprehensive_summary(const ComprehensivePerformanceReport &report)
{
    std::cout << "\n=== 综合性能分析摘要 ===" << std::endl;
    std::cout << "主要性能瓶颈: " << report.degradation_analysis.primary_bottleneck << std::endl;
    std::cout << "\n关键指标:" << std::endl;
    std::cout << "  最大线程创建开销: " << std::fixed << std::setprecision(2) << report.max_thread_creation_overhead << "x" << std::endl;
    std::cout << "  最大同步开销占比: " << (report.max_sync_overhead_ratio * 100) << "%" << std::endl;
    std::cout << "  最大False Sharing开销: " << report.max_false_sharing_overhead << "x" << std::endl;
    std::cout << "  峰值内存带宽: " << report.memory_analysis.peak_bandwidth_gbps << " GB/s" << std::endl;

    std::cout << "\n优化建议:" << std::endl;
    for (size_t i = 0; i < report.degradation_analysis.recommendations.size(); ++i)
    {
        std::cout << "  " << (i + 1) << ". " << report.degradation_analysis.recommendations[i] << std::endl;
    }

    std::cout << "\n=== 分析完成 ===" << std::endl;
}

/**
 * @brief 生成CSV格式的详细报告
 */
void export_detailed_csv_report(const ComprehensivePerformanceReport &report, const std::string &filename)
{
    std::ofstream file(filename);

    // 线程开销数据
    file << "Thread Overhead Analysis\n";
    file << "Threads,ThreadPool_ms,ThreadCreation_ms,Overhead_Ratio\n";
    for (const auto &result : report.thread_overhead_results)
    {
        file << result.num_threads << "," << result.thread_pool_time_ms << ","
             << result.thread_creation_time_ms << "," << result.overhead_ratio << "\n";
    }

    // 同步开销数据
    file << "\nSync Overhead Analysis\n";
    file << "Threads,Total_ms,Computation_ms,Sync_ms,Sync_Ratio,Avg_Barrier_ms\n";
    for (const auto &result : report.sync_overhead_results)
    {
        file << result.num_threads << "," << result.total_time_ms << "," << result.computation_time_ms
             << "," << result.sync_time_ms << "," << result.sync_overhead_ratio << "," << result.avg_barrier_time_ms << "\n";
    }

    // 内存带宽数据
    file << "\nMemory Bandwidth Analysis\n";
    file << "Threads,Bandwidth_GBps,Latency_ns,Time_ms,Efficiency_Ratio\n";
    for (const auto &result : report.memory_analysis.sequential_results)
    {
        file << result.num_threads << "," << result.bandwidth_gbps << "," << result.latency_ns
             << "," << result.time_ms << "," << result.efficiency_ratio << "\n";
    }

    file.close();
    std::cout << "详细报告已导出到: " << filename << std::endl;
}