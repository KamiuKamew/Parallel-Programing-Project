#!/usr/bin/env python3
"""
基于真实实验数据的模乘算法性能分析图表生成

数据来源:
- three_modmul_results.csv (n=4的小规模测试)
- final_modmul_results.csv (n=1024,4096,16384,32768的测试)
- modmul_performance_results.csv (n=4,1024,16384,65536,131072的测试)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# 设置matplotlib参数
plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 10
plt.rcParams["axes.titlesize"] = 12
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["xtick.labelsize"] = 9
plt.rcParams["ytick.labelsize"] = 9
plt.rcParams["legend.fontsize"] = 9


def load_real_data():
    """加载所有真实实验数据"""
    data_files = [
        "ntt/three_modmul_results.csv",
        "ntt/final_modmul_results.csv",
        "ntt/modmul_performance_results.csv",
    ]

    all_data = []

    for file_path in data_files:
        if Path(file_path).exists():
            try:
                df = pd.read_csv(file_path)
                all_data.append(df)
                print(f"Loaded {file_path}: {len(df)} records")
            except Exception as e:
                print(f"Error loading {file_path}: {e}")

    return all_data


def process_performance_data(df):
    """处理modmul_performance_results.csv数据"""
    # 提取CPU和GPU时间
    processed_data = []

    for _, row in df.iterrows():
        n = row["n"]
        algorithm = row["algorithm"]
        cpu_time = row["cpu_time_us"]
        gpu_time = row["gpu_time_us"]
        speedup = row["speedup"]
        correctness = row["correctness"]

        processed_data.append(
            {
                "n": n,
                "algorithm": algorithm,
                "cpu_time": cpu_time,
                "gpu_time": gpu_time,
                "speedup": speedup,
                "correctness": correctness,
            }
        )

    return pd.DataFrame(processed_data)


def process_final_data(df):
    """处理final_modmul_results.csv数据"""
    processed_data = []

    for _, row in df.iterrows():
        algorithm = row["Algorithm"]
        platform = row["Platform"]
        n = row["N"]
        time_us = row["Time_us"]
        correctness = row["Correctness"]

        processed_data.append(
            {
                "n": n,
                "algorithm": algorithm,
                "platform": platform,
                "time_us": time_us,
                "correctness": correctness,
            }
        )

    return pd.DataFrame(processed_data)


def create_real_performance_comparison():
    """基于真实数据创建性能对比图表"""

    # 加载数据
    data_list = load_real_data()
    if not data_list:
        print("No valid data files found!")
        return

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(
        "GPU Modular Multiplication Performance Analysis (Real Experimental Data)",
        fontsize=14,
        fontweight="bold",
    )

    # 处理modmul_performance_results.csv数据
    if len(data_list) >= 3:
        perf_data = process_performance_data(data_list[2])

        # 图1: 执行时间对比
        ax1 = axes[0, 0]

        algorithms = ["Montgomery", "Naive", "Barrett"]
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

        # 按算法分组绘制
        for i, algo in enumerate(algorithms):
            algo_data = perf_data[perf_data["algorithm"] == algo]
            if not algo_data.empty:
                sizes = algo_data["n"].values
                gpu_times = algo_data["gpu_time"].values
                ax1.plot(
                    sizes,
                    gpu_times,
                    "o-",
                    label=f"{algo}",
                    color=colors[i],
                    linewidth=2,
                    markersize=6,
                )

        ax1.set_xlabel("Problem Size (n)")
        ax1.set_ylabel("GPU Execution Time (μs)")
        ax1.set_title("GPU Algorithm Execution Time Comparison")
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # 图2: 加速比分析
        ax2 = axes[0, 1]

        for i, algo in enumerate(algorithms):
            algo_data = perf_data[perf_data["algorithm"] == algo]
            if not algo_data.empty:
                sizes = algo_data["n"].values
                speedups = algo_data["speedup"].values
                ax2.plot(
                    sizes,
                    speedups,
                    "s-",
                    label=f"{algo}",
                    color=colors[i],
                    linewidth=2,
                    markersize=6,
                )

        ax2.axhline(y=1, color="red", linestyle="--", alpha=0.7, label="CPU = GPU")
        ax2.set_xlabel("Problem Size (n)")
        ax2.set_ylabel("Speedup (CPU/GPU)")
        ax2.set_title("GPU vs CPU Speedup Analysis")
        ax2.set_xscale("log")
        ax2.grid(True, alpha=0.3)
        ax2.legend()

    # 处理final_modmul_results.csv数据
    if len(data_list) >= 2:
        final_data = process_final_data(data_list[1])

        # 图3: CPU vs GPU性能对比
        ax3 = axes[1, 0]

        # 分别处理CPU和GPU数据
        cpu_data = final_data[final_data["platform"] == "CPU"]
        gpu_data = final_data[final_data["platform"] == "GPU"]

        # CPU数据
        if not cpu_data.empty:
            cpu_naive = cpu_data[cpu_data["algorithm"] == "Naive"]
            cpu_ntt = cpu_data[cpu_data["algorithm"] == "NTT"]

            if not cpu_naive.empty:
                ax3.plot(
                    cpu_naive["n"],
                    cpu_naive["time_us"],
                    "ro-",
                    label="CPU Naive",
                    linewidth=2,
                    markersize=6,
                )
            if not cpu_ntt.empty:
                ax3.plot(
                    cpu_ntt["n"],
                    cpu_ntt["time_us"],
                    "bo-",
                    label="CPU NTT",
                    linewidth=2,
                    markersize=6,
                )

        # GPU数据
        if not gpu_data.empty:
            gpu_algorithms = gpu_data["algorithm"].unique()
            gpu_colors = ["#ff7f0e", "#2ca02c", "#d62728"]

            for i, algo in enumerate(gpu_algorithms):
                algo_data = gpu_data[gpu_data["algorithm"] == algo]
                if not algo_data.empty:
                    ax3.plot(
                        algo_data["n"],
                        algo_data["time_us"],
                        "^-",
                        label=f"GPU {algo}",
                        color=gpu_colors[i % len(gpu_colors)],
                        linewidth=2,
                        markersize=6,
                    )

        ax3.set_xlabel("Problem Size (n)")
        ax3.set_ylabel("Execution Time (μs)")
        ax3.set_title("CPU vs GPU Algorithm Performance")
        ax3.set_xscale("log")
        ax3.set_yscale("log")
        ax3.grid(True, alpha=0.3)
        ax3.legend()

    # 图4: 算法效率对比分析
    ax4 = axes[1, 1]

    # 基于small_results创建算法效率对比
    if len(data_list) >= 1:
        small_data = data_list[0]  # three_modmul_results.csv

        # 提取GPU算法数据
        gpu_rows = small_data[small_data["Platform"] == "GPU"]

        if not gpu_rows.empty:
            algorithms = gpu_rows["Algorithm"].values
            times = gpu_rows["Time_us"].values

            bars = ax4.bar(
                algorithms, times, color=["#1f77b4", "#ff7f0e", "#2ca02c"], alpha=0.7
            )
            ax4.set_ylabel("Execution Time (μs)")
            ax4.set_title("GPU Algorithm Comparison (n=4)")
            ax4.set_yscale("log")

            # 添加数值标签
            for bar, time in zip(bars, times):
                height = bar.get_height()
                ax4.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    height,
                    f"{time:.0f}μs",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )

    plt.tight_layout()
    plt.savefig("real_modmul_performance_analysis.png", dpi=300, bbox_inches="tight")
    print("Real data performance chart saved as 'real_modmul_performance_analysis.png'")
    plt.close()


def create_speedup_analysis():
    """创建加速比分析图表"""

    data_list = load_real_data()
    if len(data_list) < 3:
        print("Insufficient data for speedup analysis")
        return

    perf_data = process_performance_data(data_list[2])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle(
        "GPU Modular Multiplication Speedup Analysis (Real Data)",
        fontsize=14,
        fontweight="bold",
    )

    # 算法间加速比对比
    algorithms = ["Montgomery", "Naive", "Barrett"]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

    # 以Montgomery为基准计算相对加速比
    mont_data = perf_data[perf_data["algorithm"] == "Montgomery"]

    ax1_data = {}
    for algo in algorithms:
        algo_data = perf_data[perf_data["algorithm"] == algo]
        if not algo_data.empty:
            ax1_data[algo] = algo_data

    if "Montgomery" in ax1_data:
        for i, algo in enumerate(algorithms):
            if algo in ax1_data and algo != "Montgomery":
                algo_data = ax1_data[algo]
                mont_subset = mont_data[mont_data["n"].isin(algo_data["n"])]

                if not mont_subset.empty:
                    # 计算相对于Montgomery的加速比
                    relative_speedup = []
                    sizes = []

                    for _, row in algo_data.iterrows():
                        n = row["n"]
                        mont_row = mont_subset[mont_subset["n"] == n]
                        if not mont_row.empty:
                            mont_time = mont_row.iloc[0]["gpu_time"]
                            algo_time = row["gpu_time"]
                            if algo_time > 0:
                                speedup = mont_time / algo_time
                                relative_speedup.append(speedup)
                                sizes.append(n)

                    if relative_speedup:
                        ax1.plot(
                            sizes,
                            relative_speedup,
                            "o-",
                            label=f"{algo} vs Montgomery",
                            color=colors[i],
                            linewidth=2,
                            markersize=6,
                        )

    ax1.axhline(y=1, color="red", linestyle="--", alpha=0.7, label="Same Performance")
    ax1.set_xlabel("Problem Size (n)")
    ax1.set_ylabel("Relative Speedup")
    ax1.set_title("GPU Algorithm Relative Performance")
    ax1.set_xscale("log")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # CPU vs GPU加速比
    for i, algo in enumerate(algorithms):
        algo_data = perf_data[perf_data["algorithm"] == algo]
        if not algo_data.empty:
            valid_data = algo_data[algo_data["speedup"] > 0]
            if not valid_data.empty:
                ax2.plot(
                    valid_data["n"],
                    valid_data["speedup"],
                    "o-",
                    label=f"{algo}",
                    color=colors[i],
                    linewidth=2,
                    markersize=6,
                )

    ax2.axhline(y=1, color="red", linestyle="--", alpha=0.7, label="CPU = GPU")
    ax2.set_xlabel("Problem Size (n)")
    ax2.set_ylabel("Speedup (CPU/GPU)")
    ax2.set_title("CPU vs GPU Speedup by Algorithm")
    ax2.set_xscale("log")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("real_speedup_analysis.png", dpi=300, bbox_inches="tight")
    print("Real speedup analysis chart saved as 'real_speedup_analysis.png'")
    plt.close()


def print_data_summary():
    """打印真实数据摘要"""
    print("\n" + "=" * 50)
    print("REAL EXPERIMENTAL DATA SUMMARY")
    print("=" * 50)

    data_list = load_real_data()

    for i, df in enumerate(data_list):
        print(f"\nDataset {i+1}:")
        print(f"Shape: {df.shape}")
        print(f"Columns: {list(df.columns)}")

        if "n" in df.columns:
            sizes = sorted(df["n"].unique()) if "n" in df.columns else []
            print(f"Problem sizes: {sizes}")
        elif "N" in df.columns:
            sizes = sorted(df["N"].unique())
            print(f"Problem sizes: {sizes}")

        if "Algorithm" in df.columns:
            algos = df["Algorithm"].unique()
            print(f"Algorithms: {list(algos)}")
        elif "algorithm" in df.columns:
            algos = df["algorithm"].unique()
            print(f"Algorithms: {list(algos)}")

    print("\n" + "=" * 50)


if __name__ == "__main__":
    print("Generating charts based on REAL experimental data only...")

    # 打印数据摘要
    print_data_summary()

    # 生成图表
    create_real_performance_comparison()
    create_speedup_analysis()

    print("\nAll charts generated based on real experimental data!")
    print("Generated files:")
    print("- real_modmul_performance_analysis.png")
    print("- real_speedup_analysis.png")
