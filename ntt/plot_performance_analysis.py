#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
混合并行性能退化分析可视化脚本
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import seaborn as sns

# 设置中文字体
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "SimHei", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

# 设置绘图风格
plt.style.use("seaborn-v0_8")
sns.set_palette("husl")


def load_and_analyze_data():
    """加载并分析性能数据"""
    # 从之前的测试结果创建数据
    data = {
        "threads": [1, 2, 3, 4, 5, 6, 7, 8],
        "avg_time_us": [209845, 151334, 129828, 123242, 121970, 120977, 120709, 109871],
        "speedup": [1.0, 1.38663, 1.61633, 1.70271, 1.72047, 1.73459, 1.73843, 1.90992],
        "efficiency": [
            100.0,
            69.3316,
            53.8778,
            42.5678,
            34.4094,
            28.9098,
            24.8347,
            23.8739,
        ],
    }

    # 添加理论理想加速比和效率
    data["ideal_speedup"] = data["threads"]
    data["ideal_efficiency"] = [100.0] * len(data["threads"])

    # 添加从之前测试得到的数据：线程开销和内存带宽
    thread_overhead_data = {
        "threads": [1, 2, 4, 8],
        "overhead_time_us": [188.353, 99.39, 230.591, 261.739],
    }

    memory_bandwidth_data = {
        "threads": [1, 2, 4, 8],
        "bandwidth_gb_s": [14.8583, 22.6607, 9.57884, 8.10997],
    }

    return (
        pd.DataFrame(data),
        pd.DataFrame(thread_overhead_data),
        pd.DataFrame(memory_bandwidth_data),
    )


def analyze_task_granularity():
    """分析任务粒度问题"""
    # NTT各轮次的任务数量（从测试输出中获得）
    rounds = list(range(1, 18))
    mid_values = [2**i for i in range(17)]
    num_tasks = [131072 // (2 * mid) for mid in mid_values]
    work_per_task = mid_values

    return pd.DataFrame(
        {
            "round": rounds,
            "mid": mid_values,
            "num_tasks": num_tasks,
            "work_per_task": work_per_task,
        }
    )


def create_comprehensive_analysis():
    """创建综合性能分析图"""
    df, thread_df, memory_df = load_and_analyze_data()
    task_df = analyze_task_granularity()

    # 创建2x3的子图布局
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle("混合并行性能退化综合分析", fontsize=16, fontweight="bold")

    # 1. 性能对比图
    ax1 = axes[0, 0]
    ax1.plot(
        df["threads"],
        df["avg_time_us"] / 1000,
        "o-",
        linewidth=2,
        markersize=6,
        label="实际执行时间",
    )
    ax1.set_xlabel("线程数")
    ax1.set_ylabel("执行时间 (ms)")
    ax1.set_title("执行时间随线程数变化")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 2. 加速比对比图
    ax2 = axes[0, 1]
    ax2.plot(
        df["threads"],
        df["speedup"],
        "o-",
        linewidth=2,
        markersize=6,
        label="实际加速比",
    )
    ax2.plot(
        df["threads"],
        df["ideal_speedup"],
        "--",
        linewidth=2,
        alpha=0.7,
        label="理想加速比",
    )
    ax2.set_xlabel("线程数")
    ax2.set_ylabel("加速比")
    ax2.set_title("加速比对比")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # 3. 并行效率图
    ax3 = axes[0, 2]
    ax3.plot(
        df["threads"],
        df["efficiency"],
        "o-",
        linewidth=2,
        markersize=6,
        label="实际效率",
    )
    ax3.axhline(y=100, color="red", linestyle="--", alpha=0.7, label="理想效率")
    ax3.axhline(y=80, color="orange", linestyle=":", alpha=0.7, label="良好效率阈值")
    ax3.set_xlabel("线程数")
    ax3.set_ylabel("并行效率 (%)")
    ax3.set_title("并行效率随线程数变化")
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # 4. 任务粒度分析
    ax4 = axes[1, 0]
    # 只显示前15轮的数据，避免图表过于拥挤
    relevant_rounds = task_df[task_df["num_tasks"] >= 1].head(15)

    ax4_twin = ax4.twinx()
    bars1 = ax4.bar(
        relevant_rounds["round"],
        relevant_rounds["num_tasks"],
        alpha=0.7,
        label="可并行任务数",
        color="skyblue",
    )
    line1 = ax4_twin.plot(
        relevant_rounds["round"],
        relevant_rounds["work_per_task"],
        "ro-",
        label="每任务工作量",
    )

    ax4.set_xlabel("NTT轮次")
    ax4.set_ylabel("可并行任务数", color="blue")
    ax4_twin.set_ylabel("每任务工作量", color="red")
    ax4.set_title("NTT各轮次任务粒度分析")
    ax4.tick_params(axis="y", labelcolor="blue")
    ax4_twin.tick_params(axis="y", labelcolor="red")
    ax4.grid(True, alpha=0.3)

    # 5. 线程开销分析
    ax5 = axes[1, 1]
    ax5.plot(
        thread_df["threads"],
        thread_df["overhead_time_us"],
        "o-",
        linewidth=2,
        markersize=6,
        color="red",
        label="线程创建开销",
    )
    ax5.set_xlabel("线程数")
    ax5.set_ylabel("开销时间 (μs)")
    ax5.set_title("线程创建开销分析")
    ax5.grid(True, alpha=0.3)
    ax5.legend()

    # 6. 内存带宽竞争分析
    ax6 = axes[1, 2]
    ax6.plot(
        memory_df["threads"],
        memory_df["bandwidth_gb_s"],
        "o-",
        linewidth=2,
        markersize=6,
        color="green",
        label="内存带宽",
    )
    ax6.set_xlabel("线程数")
    ax6.set_ylabel("内存带宽 (GB/s)")
    ax6.set_title("内存带宽竞争分析")
    ax6.grid(True, alpha=0.3)
    ax6.legend()

    plt.tight_layout()
    plt.savefig("comprehensive_performance_analysis.png", dpi=300, bbox_inches="tight")
    plt.show()


def create_efficiency_breakdown():
    """创建效率分解图"""
    df, _, _ = load_and_analyze_data()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("并行效率分解分析", fontsize=14, fontweight="bold")

    # 效率随线程数变化的散点图
    ax1.scatter(df["threads"], df["efficiency"], s=100, alpha=0.7, color="red")
    ax1.plot(df["threads"], df["efficiency"], "--", alpha=0.5, color="red")

    # 标注关键点
    for i, (t, e) in enumerate(zip(df["threads"], df["efficiency"])):
        if t in [1, 2, 4, 8]:
            ax1.annotate(
                f"{e:.1f}%",
                (t, e),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=9,
            )

    ax1.axhline(y=80, color="orange", linestyle=":", alpha=0.7, label="良好效率阈值")
    ax1.set_xlabel("线程数")
    ax1.set_ylabel("并行效率 (%)")
    ax1.set_title("并行效率随线程数变化")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 效率损失分析
    efficiency_loss = [100 - e for e in df["efficiency"]]
    ax2.bar(df["threads"], efficiency_loss, alpha=0.7, color="orange")
    ax2.set_xlabel("线程数")
    ax2.set_ylabel("效率损失 (%)")
    ax2.set_title("并行效率损失分析")
    ax2.grid(True, alpha=0.3)

    # 标注损失值
    for i, (t, loss) in enumerate(zip(df["threads"], efficiency_loss)):
        if loss > 5:  # 只标注显著的损失
            ax2.text(t, loss + 1, f"{loss:.1f}%", ha="center", fontsize=9)

    plt.tight_layout()
    plt.savefig("efficiency_breakdown.png", dpi=300, bbox_inches="tight")
    plt.show()


def create_performance_heatmap():
    """创建性能热力图"""
    # 创建更详细的数据用于热力图
    threads_range = np.arange(1, 9)

    # 基于实际数据创建热力图数据
    performance_matrix = []
    efficiency_data = [100.0, 69.33, 53.88, 42.57, 34.41, 28.91, 24.83, 23.87]

    # 创建一个简化的热力图
    data_matrix = np.array(efficiency_data).reshape(1, -1)

    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(data_matrix, cmap="RdYlGn", aspect="auto", vmin=0, vmax=100)

    # 设置标签
    ax.set_xticks(range(len(threads_range)))
    ax.set_xticklabels(threads_range)
    ax.set_yticks([0])
    ax.set_yticklabels(["OpenMP NTT"])
    ax.set_xlabel("线程数")
    ax.set_title("OpenMP线程数性能效率热力图")

    # 添加文本标注
    for i, efficiency in enumerate(efficiency_data):
        ax.text(
            i,
            0,
            f"{efficiency:.1f}%",
            ha="center",
            va="center",
            color="white" if efficiency < 50 else "black",
            fontweight="bold",
        )

    # 添加颜色条
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("并行效率 (%)")

    plt.tight_layout()
    plt.savefig("performance_heatmap.png", dpi=300, bbox_inches="tight")
    plt.show()


def main():
    """主函数"""
    print("正在生成混合并行性能退化分析图表...")

    # 生成各种分析图
    create_comprehensive_analysis()
    create_efficiency_breakdown()
    create_performance_heatmap()

    print("图表生成完成！")
    print("生成的文件：")
    print("- comprehensive_performance_analysis.png")
    print("- efficiency_breakdown.png")
    print("- performance_heatmap.png")


if __name__ == "__main__":
    main()
