#!/usr/bin/env python3
"""
性能分析数据可视化脚本
生成线程数 vs 加速比和开销占比的分析图表
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pathlib import Path
import argparse

# 设置中文字体和样式
plt.rcParams["font.sans-serif"] = [
    "SimHei",
    "DejaVu Sans",
    "Arial Unicode MS",
    "sans-serif",
]
plt.rcParams["axes.unicode_minus"] = False
sns.set_style("whitegrid")
plt.rcParams["figure.figsize"] = (12, 8)


def parse_csv_report(csv_file):
    """解析CSV报告文件"""
    data = {"thread_overhead": [], "sync_overhead": [], "memory_bandwidth": []}

    with open(csv_file, "r") as f:
        lines = f.readlines()

    current_section = None
    for line in lines:
        line = line.strip()
        if not line:
            continue

        if "Thread Overhead Analysis" in line:
            current_section = "thread_overhead"
            continue
        elif "Sync Overhead Analysis" in line:
            current_section = "sync_overhead"
            continue
        elif "Memory Bandwidth Analysis" in line:
            current_section = "memory_bandwidth"
            continue
        elif line.startswith("Threads,"):
            continue  # Skip header line
        elif current_section and "," in line:
            parts = line.split(",")
            if current_section == "thread_overhead":
                data[current_section].append(
                    {
                        "threads": int(parts[0]),
                        "thread_pool_ms": float(parts[1]),
                        "thread_creation_ms": float(parts[2]),
                        "overhead_ratio": float(parts[3]),
                    }
                )
            elif current_section == "sync_overhead":
                data[current_section].append(
                    {
                        "threads": int(parts[0]),
                        "total_ms": float(parts[1]),
                        "computation_ms": float(parts[2]),
                        "sync_ms": float(parts[3]),
                        "sync_ratio": float(parts[4]),
                        "avg_barrier_ms": float(parts[5]),
                    }
                )
            elif current_section == "memory_bandwidth":
                data[current_section].append(
                    {
                        "threads": int(parts[0]),
                        "bandwidth_gbps": float(parts[1]),
                        "latency_ns": float(parts[2]),
                        "time_ms": float(parts[3]),
                        "efficiency_ratio": float(parts[4]),
                    }
                )

    return data


def plot_comprehensive_analysis(data, output_dir="./"):
    """生成综合性能分析图表"""

    # 创建2x2的子图布局
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle("并行性能瓶颈深入分析", fontsize=16, fontweight="bold")

    # 1. 线程创建开销分析
    if data["thread_overhead"]:
        threads = [d["threads"] for d in data["thread_overhead"]]
        pool_times = [d["thread_pool_ms"] for d in data["thread_overhead"]]
        creation_times = [d["thread_creation_ms"] for d in data["thread_overhead"]]
        overhead_ratios = [d["overhead_ratio"] for d in data["thread_overhead"]]

        x = np.arange(len(threads))
        width = 0.35

        bars1 = ax1.bar(
            x - width / 2,
            pool_times,
            width,
            label="线程池复用",
            alpha=0.8,
            color="skyblue",
        )
        bars2 = ax1.bar(
            x + width / 2,
            creation_times,
            width,
            label="重复创建",
            alpha=0.8,
            color="salmon",
        )

        ax1.set_xlabel("线程数")
        ax1.set_ylabel("执行时间 (ms)")
        ax1.set_title("线程创建开销对比")
        ax1.set_xticks(x)
        ax1.set_xticklabels(threads)
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # 在柱状图上标注开销比例
        for i, ratio in enumerate(overhead_ratios):
            ax1.text(
                i,
                max(creation_times[i], pool_times[i]) * 1.05,
                f"{ratio:.2f}x",
                ha="center",
                va="bottom",
                fontweight="bold",
            )

    # 2. 同步开销占比分析
    if data["sync_overhead"]:
        threads = [d["threads"] for d in data["sync_overhead"]]
        sync_ratios = [d["sync_ratio"] * 100 for d in data["sync_overhead"]]
        computation_ratios = [
            (1 - d["sync_ratio"]) * 100 for d in data["sync_overhead"]
        ]

        ax2.plot(
            threads,
            sync_ratios,
            "o-",
            linewidth=2,
            markersize=8,
            label="同步开销占比",
            color="red",
        )
        ax2.plot(
            threads,
            computation_ratios,
            "s-",
            linewidth=2,
            markersize=8,
            label="计算时间占比",
            color="blue",
        )

        ax2.set_xlabel("线程数")
        ax2.set_ylabel("时间占比 (%)")
        ax2.set_title("同步开销随线程数变化")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # 添加趋势线
        if len(threads) > 2:
            z = np.polyfit(threads, sync_ratios, 1)
            p = np.poly1d(z)
            ax2.plot(threads, p(threads), "--", alpha=0.5, color="red")

    # 3. 内存带宽竞争效率
    if data["memory_bandwidth"]:
        threads = [d["threads"] for d in data["memory_bandwidth"]]
        efficiency = [d["efficiency_ratio"] * 100 for d in data["memory_bandwidth"]]
        bandwidth = [d["bandwidth_gbps"] for d in data["memory_bandwidth"]]

        # 双y轴图
        ax3_twin = ax3.twinx()

        line1 = ax3.plot(
            threads,
            efficiency,
            "o-",
            linewidth=2,
            markersize=8,
            color="green",
            label="并行效率",
        )
        line2 = ax3_twin.plot(
            threads,
            bandwidth,
            "s-",
            linewidth=2,
            markersize=8,
            color="orange",
            label="内存带宽",
        )

        ax3.set_xlabel("线程数")
        ax3.set_ylabel("并行效率 (%)", color="green")
        ax3_twin.set_ylabel("内存带宽 (GB/s)", color="orange")
        ax3.set_title("内存带宽竞争分析")

        # 理想效率线
        ideal_efficiency = [100] * len(threads)
        ax3.plot(
            threads, ideal_efficiency, "--", alpha=0.5, color="gray", label="理想效率"
        )

        # 合并图例
        lines1, labels1 = ax3.get_legend_handles_labels()
        lines2, labels2 = ax3_twin.get_legend_handles_labels()
        ax3.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

        ax3.grid(True, alpha=0.3)

    # 4. 综合性能瓶颈分析
    # 创建一个雷达图显示不同瓶颈的严重程度
    categories = []
    values = []

    if data["thread_overhead"]:
        max_overhead = max([d["overhead_ratio"] for d in data["thread_overhead"]])
        thread_severity = min(100, (max_overhead - 1) * 50)  # 转换为0-100分数
        categories.append("线程创建开销")
        values.append(thread_severity)

    if data["sync_overhead"]:
        max_sync_ratio = max([d["sync_ratio"] for d in data["sync_overhead"]])
        sync_severity = min(100, max_sync_ratio * 500)  # 转换为0-100分数
        categories.append("同步开销")
        values.append(sync_severity)

    if data["memory_bandwidth"]:
        min_efficiency = min(
            [
                d["efficiency_ratio"]
                for d in data["memory_bandwidth"]
                if d["threads"] > 1
            ]
        )
        memory_severity = (1 - min_efficiency) * 100
        categories.append("内存带宽竞争")
        values.append(memory_severity)

    if categories and values:
        # 简化的条形图而不是雷达图
        bars = ax4.bar(
            categories, values, color=["lightcoral", "lightblue", "lightgreen"]
        )
        ax4.set_ylabel("瓶颈严重程度 (0-100)")
        ax4.set_title("性能瓶颈综合评估")
        ax4.set_ylim(0, 100)

        # 添加数值标签
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax4.text(
                bar.get_x() + bar.get_width() / 2.0,
                height + 1,
                f"{value:.1f}",
                ha="center",
                va="bottom",
                fontweight="bold",
            )

        # 添加严重程度区间指示
        ax4.axhline(y=30, color="green", linestyle="--", alpha=0.7, label="轻微瓶颈")
        ax4.axhline(y=60, color="orange", linestyle="--", alpha=0.7, label="中等瓶颈")
        ax4.axhline(y=80, color="red", linestyle="--", alpha=0.7, label="严重瓶颈")
        ax4.legend()

    plt.tight_layout()

    # 保存图表
    output_path = Path(output_dir) / "comprehensive_performance_analysis.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"综合性能分析图表已保存到: {output_path}")

    return fig


def generate_performance_summary_table(data, output_dir="./"):
    """生成性能分析摘要表格"""

    summary_data = []

    # 线程创建开销摘要
    if data["thread_overhead"]:
        max_overhead = max([d["overhead_ratio"] for d in data["thread_overhead"]])
        summary_data.append(
            {
                "性能指标": "线程创建开销",
                "最大值": f"{max_overhead:.2f}x",
                "建议阈值": "< 1.5x",
                "状态": "⚠️ 需要优化" if max_overhead > 1.5 else "✅ 正常",
            }
        )

    # 同步开销摘要
    if data["sync_overhead"]:
        max_sync = max([d["sync_ratio"] for d in data["sync_overhead"]])
        summary_data.append(
            {
                "性能指标": "同步开销占比",
                "最大值": f"{max_sync*100:.1f}%",
                "建议阈值": "< 20%",
                "状态": "⚠️ 需要优化" if max_sync > 0.2 else "✅ 正常",
            }
        )

    # 内存效率摘要
    if data["memory_bandwidth"]:
        min_efficiency = min(
            [
                d["efficiency_ratio"]
                for d in data["memory_bandwidth"]
                if d["threads"] > 1
            ]
        )
        summary_data.append(
            {
                "性能指标": "最低并行效率",
                "最大值": f"{min_efficiency*100:.1f}%",
                "建议阈值": "> 80%",
                "状态": "⚠️ 需要优化" if min_efficiency < 0.8 else "✅ 正常",
            }
        )

    # 创建DataFrame并保存
    df = pd.DataFrame(summary_data)

    # 生成表格图片
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis("tight")
    ax.axis("off")

    table = ax.table(
        cellText=df.values.tolist(),
        colLabels=df.columns.tolist(),
        cellLoc="center",
        loc="center",
        colWidths=[0.3, 0.2, 0.2, 0.3],
    )

    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.5)

    # 设置表格样式
    for i in range(len(df.columns)):
        table[(0, i)].set_facecolor("#40466e")
        table[(0, i)].set_text_props(weight="bold", color="white")

    for i in range(1, len(df) + 1):
        for j in range(len(df.columns)):
            if "⚠️" in str(df.iloc[i - 1, j]):
                table[(i, j)].set_facecolor("#ffcccc")
            elif "✅" in str(df.iloc[i - 1, j]):
                table[(i, j)].set_facecolor("#ccffcc")

    plt.title("性能分析摘要表", fontsize=16, fontweight="bold", pad=20)

    # 保存表格
    output_path = Path(output_dir) / "performance_summary_table.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"性能摘要表已保存到: {output_path}")

    return df


def main():
    parser = argparse.ArgumentParser(description="性能分析数据可视化")
    parser.add_argument(
        "--csv",
        type=str,
        default="performance_analysis_test.csv",
        help="CSV报告文件路径",
    )
    parser.add_argument("--output", type=str, default="./", help="输出目录")

    args = parser.parse_args()

    if not Path(args.csv).exists():
        print(f"错误: CSV文件 {args.csv} 不存在")
        return

    print(f"解析CSV报告: {args.csv}")
    data = parse_csv_report(args.csv)

    print("生成综合性能分析图表...")
    plot_comprehensive_analysis(data, args.output)

    print("生成性能摘要表...")
    summary_df = generate_performance_summary_table(data, args.output)

    print("\n=== 性能分析摘要 ===")
    print(summary_df.to_string(index=False))

    print(f"\n所有图表已保存到目录: {args.output}")


if __name__ == "__main__":
    main()
