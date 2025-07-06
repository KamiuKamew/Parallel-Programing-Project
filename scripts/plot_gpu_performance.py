#!/usr/bin/env python3
"""
GPU NTT优化效果性能分析脚本
生成性能对比图表和详细分析报告
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os

# 设置中文字体支持
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "SimHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False


def create_sample_data():
    """创建示例性能数据（如果没有实际测试数据）"""
    sizes = [1024, 4096, 16384, 65536, 131072, 262144, 524288]

    # 模拟性能数据（基于Lab5.tex中的分析）
    cpu_times = [50, 200, 1000, 8500, 37000, 85000, 190000]

    gpu_basic_times = [2000, 2500, 3000, 11000, 35000, 63000, 137000]

    gpu_optimized_times = [1800, 2000, 2200, 7000, 22000, 40000, 85000]

    data = []
    for i, n in enumerate(sizes):
        basic_speedup = cpu_times[i] / gpu_basic_times[i]
        opt_speedup = cpu_times[i] / gpu_optimized_times[i]
        improvement = gpu_basic_times[i] / gpu_optimized_times[i]

        data.append(
            {
                "Size": n,
                "CPU_Time_us": cpu_times[i],
                "GPU_Basic_Time_us": gpu_basic_times[i],
                "GPU_Optimized_Time_us": gpu_optimized_times[i],
                "Basic_Speedup": basic_speedup,
                "Optimized_Speedup": opt_speedup,
                "Improvement_Ratio": improvement,
            }
        )

    return pd.DataFrame(data)


def load_performance_data():
    """加载性能测试数据"""
    csv_file = "gpu_optimization_benchmark.csv"

    if os.path.exists(csv_file):
        print(f"Loading data from {csv_file}")
        return pd.read_csv(csv_file)
    else:
        print("No benchmark data found, creating sample data...")
        return create_sample_data()


def plot_performance_comparison(df):
    """绘制性能对比图"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

    # 图1: 执行时间对比（对数坐标）
    ax1.loglog(
        df["Size"], df["CPU_Time_us"], "bo-", label="CPU", linewidth=2, markersize=6
    )
    ax1.loglog(
        df["Size"],
        df["GPU_Basic_Time_us"],
        "rs-",
        label="GPU Basic",
        linewidth=2,
        markersize=6,
    )
    ax1.loglog(
        df["Size"],
        df["GPU_Optimized_Time_us"],
        "g^-",
        label="GPU Optimized",
        linewidth=2,
        markersize=6,
    )
    ax1.set_xlabel("Problem Size (n)", fontsize=12)
    ax1.set_ylabel("Execution Time (μs)", fontsize=12)
    ax1.set_title(
        "Execution Time Comparison (Log-Log Scale)", fontsize=14, fontweight="bold"
    )
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)

    # 图2: 加速比对比
    ax2.semilogx(
        df["Size"],
        df["Basic_Speedup"],
        "rs-",
        label="GPU Basic Speedup",
        linewidth=2,
        markersize=6,
    )
    ax2.semilogx(
        df["Size"],
        df["Optimized_Speedup"],
        "g^-",
        label="GPU Optimized Speedup",
        linewidth=2,
        markersize=6,
    )
    ax2.axhline(y=1, color="k", linestyle="--", alpha=0.5, label="No Speedup")
    ax2.set_xlabel("Problem Size (n)", fontsize=12)
    ax2.set_ylabel("Speedup vs CPU", fontsize=12)
    ax2.set_title("GPU Speedup vs CPU", fontsize=14, fontweight="bold")
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    # 图3: 优化提升效果
    ax3.semilogx(df["Size"], df["Improvement_Ratio"], "mo-", linewidth=2, markersize=8)
    ax3.axhline(y=1, color="k", linestyle="--", alpha=0.5, label="No Improvement")
    ax3.set_xlabel("Problem Size (n)", fontsize=12)
    ax3.set_ylabel("Improvement Ratio", fontsize=12)
    ax3.set_title(
        "GPU Optimization Improvement\n(Optimized / Basic)",
        fontsize=14,
        fontweight="bold",
    )
    ax3.legend(fontsize=11)
    ax3.grid(True, alpha=0.3)

    # 图4: 性能效率分析
    efficiency_basic = df["Basic_Speedup"] / (df["Size"] / 1000)  # 标准化效率
    efficiency_opt = df["Optimized_Speedup"] / (df["Size"] / 1000)

    ax4.semilogx(
        df["Size"],
        efficiency_basic,
        "rs-",
        label="Basic Efficiency",
        linewidth=2,
        markersize=6,
    )
    ax4.semilogx(
        df["Size"],
        efficiency_opt,
        "g^-",
        label="Optimized Efficiency",
        linewidth=2,
        markersize=6,
    )
    ax4.set_xlabel("Problem Size (n)", fontsize=12)
    ax4.set_ylabel("Efficiency (Speedup per 1K elements)", fontsize=12)
    ax4.set_title("GPU Efficiency Analysis", fontsize=14, fontweight="bold")
    ax4.legend(fontsize=11)
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("gpu_performance_analysis.png", dpi=300, bbox_inches="tight")
    print("Performance comparison chart saved: gpu_performance_analysis.png")


def generate_performance_report(df):
    """生成详细的性能分析报告"""
    report = []
    report.append("# GPU NTT优化效果分析报告")
    report.append("=" * 50)
    report.append("")

    # 基本统计
    max_basic_speedup = df["Basic_Speedup"].max()
    max_opt_speedup = df["Optimized_Speedup"].max()
    max_improvement = df["Improvement_Ratio"].max()

    max_basic_idx = df["Basic_Speedup"].idxmax()
    max_opt_idx = df["Optimized_Speedup"].idxmax()
    max_imp_idx = df["Improvement_Ratio"].idxmax()

    report.append("## 关键性能指标")
    report.append(
        f"- 最佳基础GPU加速比: {max_basic_speedup:.3f}x (n={df.loc[max_basic_idx, 'Size']})"
    )
    report.append(
        f"- 最佳优化GPU加速比: {max_opt_speedup:.3f}x (n={df.loc[max_opt_idx, 'Size']})"
    )
    report.append(
        f"- 最大优化提升倍数: {max_improvement:.3f}x (n={df.loc[max_imp_idx, 'Size']})"
    )
    report.append("")

    # 大规模问题分析
    large_scale = df[df["Size"] >= 65536]
    if not large_scale.empty:
        avg_basic_speedup = large_scale["Basic_Speedup"].mean()
        avg_opt_speedup = large_scale["Optimized_Speedup"].mean()
        avg_improvement = large_scale["Improvement_Ratio"].mean()

        report.append("## 大规模问题性能 (n≥65536)")
        report.append(f"- 平均基础GPU加速比: {avg_basic_speedup:.3f}x")
        report.append(f"- 平均优化GPU加速比: {avg_opt_speedup:.3f}x")
        report.append(f"- 平均优化提升: {avg_improvement:.3f}x")
        report.append("")

    # 详细数据表
    report.append("## 详细测试结果")
    report.append(
        "| Size | CPU(μs) | GPU基础(μs) | GPU优化(μs) | 基础加速比 | 优化加速比 | 优化提升 |"
    )
    report.append(
        "|------|---------|-------------|-------------|-----------|-----------|----------|"
    )

    for _, row in df.iterrows():
        report.append(
            f"| {row['Size']:6d} | {row['CPU_Time_us']:7.0f} | "
            f"{row['GPU_Basic_Time_us']:9.0f} | {row['GPU_Optimized_Time_us']:9.0f} | "
            f"{row['Basic_Speedup']:8.3f} | {row['Optimized_Speedup']:8.3f} | "
            f"{row['Improvement_Ratio']:7.3f} |"
        )

    report.append("")

    # 优化策略分析
    report.append("## 优化策略效果分析")

    if max_opt_speedup > max_basic_speedup:
        report.append("✅ **优化成功**: GPU优化版本显著提升了性能")
    else:
        report.append("⚠️  **需要改进**: 优化效果不明显，需要进一步调优")

    if avg_improvement > 1.2:  # 假设20%以上提升为显著
        report.append("✅ **优化策略有效**: 平均优化提升超过20%")
    else:
        report.append("📊 **优化空间**: 还有进一步优化的空间")

    # 保存报告
    with open("gpu_performance_report.md", "w", encoding="utf-8") as f:
        f.write("\n".join(report))

    print("Performance report saved: gpu_performance_report.md")

    # 打印摘要到控制台
    print("\n" + "=" * 50)
    print("GPU优化效果摘要:")
    print(f"最佳优化加速比: {max_opt_speedup:.3f}x")
    print(f"最大优化提升: {max_improvement:.3f}x")
    if not large_scale.empty:
        print(f"大规模平均提升: {avg_improvement:.3f}x")
    print("=" * 50)


def main():
    """主函数"""
    print("GPU NTT优化效果分析工具")
    print("=" * 30)

    # 加载数据
    df = load_performance_data()

    if df.empty:
        print("错误: 无法加载性能数据")
        sys.exit(1)

    print(f"加载了 {len(df)} 个测试数据点")

    # 生成图表
    plot_performance_comparison(df)

    # 生成报告
    generate_performance_report(df)

    print("\n分析完成! 生成的文件:")
    print("- gpu_performance_analysis.png: 性能对比图表")
    print("- gpu_performance_report.md: 详细分析报告")


if __name__ == "__main__":
    main()
