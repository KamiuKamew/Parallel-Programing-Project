#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings

warnings.filterwarnings("ignore")

# 设置matplotlib参数，按照规则使用英文
plt.rcParams["font.size"] = 12
plt.rcParams["figure.figsize"] = (10, 8)
plt.rcParams["axes.grid"] = True
plt.rcParams["grid.alpha"] = 0.3


def load_data():
    """加载和处理数据"""
    df = pd.read_csv("fixed_modmul_results.csv")
    df["Algorithm_Type"] = (
        df["Algorithm"].str.replace("CPU ", "").str.replace("GPU ", "")
    )

    modulus_names = {
        7340033: "7.3M",
        104857601: "104M",
        469762049: "469M",
        998244353: "998M",
    }
    df["Modulus_Name"] = df["Modulus"].map(modulus_names)
    return df


def plot_gpu_vs_cpu_speedup():
    """绘制GPU vs CPU性能分析图"""
    df = load_data()

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    cpu_data = df[df["Platform"] == "CPU"]
    gpu_mont_data = df[
        (df["Platform"] == "GPU") & (df["Algorithm_Type"] == "Montgomery")
    ]

    # 选择代表性模数
    modulus = 104857601
    cpu_subset = cpu_data[cpu_data["Modulus"] == modulus]
    gpu_subset = gpu_mont_data[gpu_mont_data["Modulus"] == modulus]

    merged = pd.merge(
        cpu_subset[["N", "Time_us"]],
        gpu_subset[["N", "Time_us"]],
        on="N",
        suffixes=("_cpu", "_gpu"),
    )
    merged["Speedup"] = merged["Time_us_cpu"] / merged["Time_us_gpu"]

    # 绘制主要曲线
    ax.plot(
        merged["N"],
        merged["Speedup"],
        "o-",
        linewidth=4,
        markersize=10,
        color="#2E86C1",
        label="GPU Montgomery vs CPU Montgomery",
    )

    # 添加关键参考线
    ax.axhline(
        y=1, color="red", linestyle="--", alpha=0.8, linewidth=3, label="Breakeven Line"
    )
    ax.axvline(
        x=16384,
        color="orange",
        linestyle="--",
        alpha=0.8,
        linewidth=3,
        label="Performance Crossover (n=16384)",
    )

    # 添加数据标签
    for i, row in merged.iterrows():
        ax.annotate(
            f'{row["Speedup"]:.2f}x',
            (row["N"], row["Speedup"]),
            textcoords="offset points",
            xytext=(0, 15),
            ha="center",
            fontweight="bold",
            fontsize=11,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )

    # 添加性能阶段文字标注
    ax.text(
        100,
        0.3,
        "Small Scale:\nGPU Overhead\nDominates",
        ha="center",
        va="center",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.8),
    )
    ax.text(
        65536,
        2.2,
        "Large Scale:\nGPU Advantage\nClear",
        ha="center",
        va="center",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.8),
    )

    ax.set_xlabel("Problem Size (n)", fontsize=14, fontweight="bold")
    ax.set_ylabel("GPU Speedup vs CPU", fontsize=14, fontweight="bold")
    ax.set_title(
        "GPU vs CPU Performance Analysis\n(Montgomery Algorithm, p=104M)",
        fontsize=16,
        fontweight="bold",
    )
    ax.set_xscale("log", base=2)
    ax.set_ylim(0, 2.5)
    ax.legend(fontsize=11, loc="upper left")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("gpu_vs_cpu_speedup.png", dpi=300, bbox_inches="tight")
    plt.show()
    print("Generated: gpu_vs_cpu_speedup.png")


def plot_gpu_algorithm_comparison():
    """绘制GPU算法性能对比图"""
    df = load_data()

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    # 选择中等模数进行详细分析
    gpu_data = df[(df["Platform"] == "GPU") & (df["Modulus"] == 104857601)]
    algorithms = ["Naive", "Montgomery", "Barrett"]
    colors = ["#E74C3C", "#F39C12", "#27AE60"]
    markers = ["o", "s", "^"]

    for i, algo in enumerate(algorithms):
        algo_data = gpu_data[gpu_data["Algorithm_Type"] == algo]
        ax.plot(
            algo_data["N"],
            algo_data["Time_us"],
            marker=markers[i],
            linestyle="-",
            label=f"GPU {algo}",
            color=colors[i],
            linewidth=3,
            markersize=8,
            markerfacecolor="white",
            markeredgewidth=2,
            markeredgecolor=colors[i],
        )

    ax.set_xlabel("Problem Size (n)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Execution Time (μs)", fontsize=14, fontweight="bold")
    ax.set_title(
        "GPU Algorithm Performance Comparison\n(p=104M)", fontsize=16, fontweight="bold"
    )
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.legend(fontsize=12, loc="upper left")
    ax.grid(True, alpha=0.3)

    # 添加性能特点标注
    ax.annotate(
        "Small problems:\nNaive performs best",
        xy=(4, 1624),
        xytext=(100, 5000),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.2"),
        fontsize=10,
        ha="center",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.8),
    )

    ax.annotate(
        "Large problems:\nMontgomery and Barrett\nshow optimization benefits",
        xy=(131072, 58991),
        xytext=(20000, 20000),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-0.2"),
        fontsize=10,
        ha="center",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.8),
    )

    plt.tight_layout()
    plt.savefig("gpu_algorithm_comparison.png", dpi=300, bbox_inches="tight")
    plt.show()
    print("Generated: gpu_algorithm_comparison.png")


def plot_algorithm_relative_performance():
    """绘制算法相对性能分析图"""
    df = load_data()

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    # 选择中等模数进行详细分析
    gpu_data = df[(df["Platform"] == "GPU") & (df["Modulus"] == 104857601)]
    naive_data = gpu_data[gpu_data["Algorithm_Type"] == "Naive"].set_index("N")[
        "Time_us"
    ]

    colors = ["#F39C12", "#27AE60"]
    markers = ["s", "^"]

    for i, algo in enumerate(["Montgomery", "Barrett"]):
        algo_data = gpu_data[gpu_data["Algorithm_Type"] == algo].set_index("N")[
            "Time_us"
        ]
        relative_perf = naive_data / algo_data  # >1表示比Naive快
        ax.plot(
            relative_perf.index,
            relative_perf.values,
            marker=markers[i],
            linestyle="-",
            label=f"{algo} vs Naive",
            color=colors[i],
            linewidth=3,
            markersize=8,
            markerfacecolor="white",
            markeredgewidth=2,
            markeredgecolor=colors[i],
        )

    ax.axhline(
        y=1,
        color="red",
        linestyle="--",
        alpha=0.8,
        linewidth=2,
        label="Same Performance as Naive",
    )

    # 移除了性能区域填充，保持图表简洁

    ax.set_xlabel("Problem Size (n)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Speedup vs Naive", fontsize=14, fontweight="bold")
    ax.set_title(
        "GPU Algorithm Relative Performance Analysis\n(p=104M)",
        fontsize=16,
        fontweight="bold",
    )
    ax.set_xscale("log", base=2)
    ax.set_ylim(0.4, 1.8)
    ax.legend(fontsize=11, loc="upper right")
    ax.grid(True, alpha=0.3)

    # 添加关键发现标注
    ax.annotate(
        "Barrett shows\nsignificant optimization\nat n=1024 and n=65536",
        xy=(1024, 1.587),
        xytext=(20000, 1.4),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-0.3"),
        fontsize=10,
        ha="center",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.8),
    )

    plt.tight_layout()
    plt.savefig("algorithm_relative_performance.png", dpi=300, bbox_inches="tight")
    plt.show()
    print("Generated: algorithm_relative_performance.png")


def main():
    """主函数"""
    print("🔧 Loading fixed experimental data...")

    print("📈 Generating GPU vs CPU speedup analysis...")
    plot_gpu_vs_cpu_speedup()

    print("📊 Generating GPU algorithm comparison...")
    plot_gpu_algorithm_comparison()

    print("🔍 Generating algorithm relative performance analysis...")
    plot_algorithm_relative_performance()

    print("\n✅ Individual Visualization Complete!")
    print("\nGenerated files:")
    print("  📄 gpu_vs_cpu_speedup.png")
    print("  📄 gpu_algorithm_comparison.png")
    print("  📄 algorithm_relative_performance.png")

    print("\n📊 Summary:")
    print("  • GPU vs CPU: Shows performance crossover at n=16384")
    print("  • Algorithm Comparison: Demonstrates optimization benefits")
    print("  • Relative Performance: Validates reduction algorithm improvements")


if __name__ == "__main__":
    main()
