#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MPI NTT性能测试数据可视化脚本
生成加速比、效率、可扩展性分析图表
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rcParams

# 设置英文字体和图表样式
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False
sns.set_style("whitegrid")

# 实验数据 (基于实际测试结果)
# 时间单位：微秒 (us)
performance_data = {
    "n=4": {
        "p=7340033": {
            "serial": 63.59,
            "1_process": 0.38,
            "2_process": 0.42,
            "2p_2t": 0.54,
            "2p_4t": 4.41,
        }
    },
    "n=131072": {
        "p=7340033": {
            "serial": 526.01,
            "1_process": 543.70,
            "2_process": 301.65,
            "2p_2t": 304.41,
            "2p_4t": 365.87,
        },
        "p=104857601": {
            "serial": 458.33,
            "1_process": 489.73,
            "2_process": 276.15,
            "2p_2t": 284.58,
            "2p_4t": 296.37,
        },
        "p=469762049": {
            "serial": 496.61,
            "1_process": 453.01,
            "2_process": 268.62,
            "2p_2t": 256.08,
            "2p_4t": 262.89,
        },
        "p=1337006139375617": {
            "serial": 407.49,
            "1_process": 469.16,
            "2_process": 262.91,
            "2p_2t": 251.56,
            "2p_4t": 318.96,
        },
    },
}


def calculate_speedup_efficiency():
    """计算加速比和效率"""
    results = []

    for problem_size, moduli in performance_data.items():
        for modulus, times in moduli.items():
            serial_time = times["serial"]

            for config, time in times.items():
                if config == "serial":
                    continue

                speedup = serial_time / time

                # 确定进程数和线程数
                processes = 1
                threads = 1
                if config == "1_process":
                    processes = 1
                    threads = 1
                elif config == "2_process":
                    processes = 2
                    threads = 1
                elif config == "2p_2t":
                    processes = 2
                    threads = 2
                elif config == "2p_4t":
                    processes = 2
                    threads = 4

                total_cores = processes * threads
                efficiency = speedup / total_cores * 100

                results.append(
                    {
                        "problem_size": problem_size,
                        "modulus": modulus,
                        "config": config,
                        "processes": processes,
                        "threads": threads,
                        "total_cores": total_cores,
                        "time": time,
                        "speedup": speedup,
                        "efficiency": efficiency,
                    }
                )

    return pd.DataFrame(results)


def plot_speedup_comparison():
    """绘制加速比对比图"""
    df = calculate_speedup_efficiency()

    # 筛选n=131072的数据
    large_scale = df[df["problem_size"] == "n=131072"]

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    moduli = large_scale["modulus"].unique()

    for i, mod in enumerate(moduli):
        row = i // 2
        col = i % 2

        mod_data = large_scale[large_scale["modulus"] == mod]

        # 准备数据
        configs = mod_data["config"]
        speedups = mod_data["speedup"]

        # 绘制条形图
        bars = axes[row, col].bar(
            range(len(configs)),
            speedups,
            color=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"],
        )

        # 添加数值标签
        for j, (bar, speedup) in enumerate(zip(bars, speedups)):
            axes[row, col].text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.02,
                f"{speedup:.2f}x",
                ha="center",
                va="bottom",
                fontsize=10,
            )

        # 添加理想加速比线
        ideal_speedup = [1, 2, 4, 8]
        axes[row, col].plot(
            range(len(configs)),
            ideal_speedup[: len(configs)],
            "r--",
            alpha=0.7,
            label="理想加速比",
        )

        axes[row, col].set_title(f"模数: {mod}", fontsize=12, fontweight="bold")
        axes[row, col].set_ylabel("加速比", fontsize=11)
        axes[row, col].set_xticks(range(len(configs)))
        axes[row, col].set_xticklabels(
            ["1进程", "2进程", "2进程2线程", "2进程4线程"], rotation=45, ha="right"
        )
        axes[row, col].grid(True, alpha=0.3)
        axes[row, col].legend()

    plt.tight_layout()
    plt.savefig("../image/speedup_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()


def plot_efficiency_analysis():
    """绘制并行效率分析图"""
    df = calculate_speedup_efficiency()
    large_scale = df[df["problem_size"] == "n=131072"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # 效率对比
    moduli = large_scale["modulus"].unique()
    configs = ["1_process", "2_process", "2p_2t", "2p_4t"]
    config_labels = ["1进程", "2进程", "2进程2线程", "2进程4线程"]

    x = np.arange(len(config_labels))
    width = 0.2

    for i, mod in enumerate(moduli):
        mod_data = large_scale[large_scale["modulus"] == mod]
        efficiencies = []

        for config in configs:
            eff = mod_data[mod_data["config"] == config]["efficiency"].iloc[0]
            efficiencies.append(eff)

        ax1.bar(x + i * width, efficiencies, width, label=f"模数{i+1}", alpha=0.8)

    ax1.set_xlabel("并行配置", fontsize=12)
    ax1.set_ylabel("并行效率 (%)", fontsize=12)
    ax1.set_title("不同模数下的并行效率对比", fontsize=14, fontweight="bold")
    ax1.set_xticks(x + width * 1.5)
    ax1.set_xticklabels(config_labels)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 平均效率趋势
    avg_efficiency = []
    for config in configs:
        avg_eff = large_scale[large_scale["config"] == config]["efficiency"].mean()
        avg_efficiency.append(avg_eff)

    ax2.plot(
        config_labels,
        avg_efficiency,
        "o-",
        linewidth=2,
        markersize=8,
        color="#2ca02c",
        label="平均效率",
    )
    ax2.axhline(y=100, color="r", linestyle="--", alpha=0.7, label="理想效率")
    ax2.axhline(y=80, color="orange", linestyle="--", alpha=0.7, label="良好效率阈值")

    for i, eff in enumerate(avg_efficiency):
        ax2.text(i, eff + 2, f"{eff:.1f}%", ha="center", va="bottom", fontsize=10)

    ax2.set_xlabel("并行配置", fontsize=12)
    ax2.set_ylabel("平均并行效率 (%)", fontsize=12)
    ax2.set_title("平均并行效率趋势", fontsize=14, fontweight="bold")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("../image/efficiency_analysis.png", dpi=300, bbox_inches="tight")
    plt.close()


def plot_scalability_analysis():
    """绘制可扩展性分析图"""
    df = calculate_speedup_efficiency()
    large_scale = df[df["problem_size"] == "n=131072"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # 强可扩展性 (固定问题规模)
    processes = [1, 2]
    avg_times = []

    for p in processes:
        if p == 1:
            avg_time = large_scale[large_scale["config"] == "1_process"]["time"].mean()
        else:
            avg_time = large_scale[large_scale["config"] == "2_process"]["time"].mean()
        avg_times.append(avg_time)

    ax1.plot(
        processes,
        avg_times,
        "o-",
        linewidth=2,
        markersize=8,
        color="#1f77b4",
        label="实际性能",
    )

    # 理想可扩展性
    ideal_times = [avg_times[0] / p for p in processes]
    ax1.plot(
        processes,
        ideal_times,
        "--",
        linewidth=2,
        color="red",
        alpha=0.7,
        label="理想可扩展性",
    )

    ax1.set_xlabel("进程数", fontsize=12)
    ax1.set_ylabel("平均执行时间 (μs)", fontsize=12)
    ax1.set_title("强可扩展性分析 (n=131072)", fontsize=14, fontweight="bold")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 混合并行效果对比
    hybrid_configs = ["2_process", "2p_2t", "2p_4t"]
    hybrid_labels = ["2进程", "2进程2线程", "2进程4线程"]

    avg_times_hybrid = []
    for config in hybrid_configs:
        avg_time = large_scale[large_scale["config"] == config]["time"].mean()
        avg_times_hybrid.append(avg_time)

    bars = ax2.bar(
        hybrid_labels,
        avg_times_hybrid,
        color=["#ff7f0e", "#2ca02c", "#d62728"],
        alpha=0.8,
    )

    for bar, time in zip(bars, avg_times_hybrid):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 5,
            f"{time:.1f}μs",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    ax2.set_ylabel("平均执行时间 (μs)", fontsize=12)
    ax2.set_title("混合并行配置对比", fontsize=14, fontweight="bold")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("../image/scalability_analysis.png", dpi=300, bbox_inches="tight")
    plt.close()


def create_performance_summary_table():
    """创建性能汇总表"""
    df = calculate_speedup_efficiency()

    # 生成汇总表
    summary = (
        df.groupby(["problem_size", "config"])
        .agg({"time": "mean", "speedup": "mean", "efficiency": "mean"})
        .round(2)
    )

    print("性能测试结果汇总表:")
    print("=" * 80)
    print(summary.to_string())
    print("=" * 80)

    # 保存为CSV
    summary.to_csv("../image/performance_summary.csv")

    return summary


def main():
    """主函数：生成所有图表"""
    print("正在生成性能分析图表...")

    # 生成图表
    plot_speedup_comparison()
    print("✓ 加速比对比图已保存: ../image/speedup_comparison.png")

    plot_efficiency_analysis()
    print("✓ 效率分析图已保存: ../image/efficiency_analysis.png")

    plot_scalability_analysis()
    print("✓ 可扩展性分析图已保存: ../image/scalability_analysis.png")

    # 生成汇总表
    create_performance_summary_table()
    print("✓ 性能汇总表已保存: ../image/performance_summary.csv")

    print("\n所有图表和数据表已生成完成！")


if __name__ == "__main__":
    main()
