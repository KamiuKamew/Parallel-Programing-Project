#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MPI NTT Performance Test Data Visualization Script (English Version)
Generate speedup, efficiency, and scalability analysis charts
Focus on different problem sizes rather than different moduli
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rcParams

# Set English font and chart style
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False
sns.set_style("whitegrid")

# Experimental data (based on actual test results)
# Time unit: microseconds (us)
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
    """Calculate speedup and efficiency metrics"""
    results = []

    for problem_size, moduli in performance_data.items():
        for modulus, times in moduli.items():
            serial_time = times["serial"]

            for config, time in times.items():
                if config == "serial":
                    continue

                speedup = serial_time / time

                # Determine process and thread counts
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
    """Plot speedup comparison focusing on different problem sizes"""
    df = calculate_speedup_efficiency()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Problem sizes to compare
    problem_sizes = ["n=4", "n=131072"]
    problem_labels = ["Small Scale (n=4)", "Large Scale (n=131,072)"]
    configs = ["1_process", "2_process", "2p_2t", "2p_4t"]
    config_labels = ["1 Process", "2 Processes", "2P+2T", "2P+4T"]

    x = np.arange(len(config_labels))
    width = 0.35

    # Left plot: Comparison between different problem sizes
    for i, (size, label) in enumerate(zip(problem_sizes, problem_labels)):
        size_data = df[df["problem_size"] == size]
        avg_speedups = []

        for config in configs:
            speedup = size_data[size_data["config"] == config]["speedup"].mean()
            avg_speedups.append(speedup)

        bars = ax1.bar(x + i * width, avg_speedups, width, label=label, alpha=0.8)

        # Add value labels
        for j, (bar, speedup) in enumerate(zip(bars, avg_speedups)):
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(avg_speedups) * 0.02,
                f"{speedup:.1f}x",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    ax1.set_xlabel("Parallel Configuration", fontsize=12)
    ax1.set_ylabel("Average Speedup", fontsize=12)
    ax1.set_title(
        "Speedup Comparison: Problem Size Impact", fontsize=14, fontweight="bold"
    )
    ax1.set_xticks(x + width / 2)
    ax1.set_xticklabels(config_labels)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Right plot: Focus on large-scale performance (n=131072)
    large_scale = df[df["problem_size"] == "n=131072"]
    avg_speedups_large = []

    for config in configs:
        speedup = large_scale[large_scale["config"] == config]["speedup"].mean()
        avg_speedups_large.append(speedup)

    bars2 = ax2.bar(
        config_labels,
        avg_speedups_large,
        color=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"],
        alpha=0.8,
    )

    # Add ideal speedup line
    ideal_speedup = [1, 2, 4, 8]
    ax2.plot(
        config_labels,
        ideal_speedup[: len(config_labels)],
        "r--",
        linewidth=2,
        alpha=0.7,
        label="Ideal Speedup",
    )

    # Add value labels
    for bar, speedup in zip(bars2, avg_speedups_large):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.03,
            f"{speedup:.2f}x",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    ax2.set_xlabel("Parallel Configuration", fontsize=12)
    ax2.set_ylabel("Average Speedup", fontsize=12)
    ax2.set_title("Large Scale Performance (n=131,072)", fontsize=14, fontweight="bold")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("./image/speedup_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()


def plot_efficiency_analysis():
    """Plot parallel efficiency analysis with emphasis on problem sizes"""
    df = calculate_speedup_efficiency()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Left plot: Efficiency comparison between problem sizes
    problem_sizes = ["n=4", "n=131072"]
    problem_labels = ["Small (n=4)", "Large (n=131,072)"]
    configs = ["1_process", "2_process", "2p_2t", "2p_4t"]
    config_labels = ["1 Process", "2 Processes", "2P+2T", "2P+4T"]

    x = np.arange(len(config_labels))
    width = 0.35

    for i, (size, label) in enumerate(zip(problem_sizes, problem_labels)):
        size_data = df[df["problem_size"] == size]
        avg_efficiencies = []

        for config in configs:
            eff = size_data[size_data["config"] == config]["efficiency"].mean()
            avg_efficiencies.append(eff)

        ax1.bar(x + i * width, avg_efficiencies, width, label=label, alpha=0.8)

    ax1.set_xlabel("Parallel Configuration", fontsize=12)
    ax1.set_ylabel("Parallel Efficiency (%)", fontsize=12)
    ax1.set_title(
        "Efficiency Comparison: Problem Size Impact", fontsize=14, fontweight="bold"
    )
    ax1.set_xticks(x + width / 2)
    ax1.set_xticklabels(config_labels)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Right plot: Efficiency trend for large scale
    large_scale = df[df["problem_size"] == "n=131072"]
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
        label="Average Efficiency",
    )
    ax2.axhline(y=100, color="r", linestyle="--", alpha=0.7, label="Ideal Efficiency")
    ax2.axhline(
        y=80,
        color="orange",
        linestyle="--",
        alpha=0.7,
        label="Good Efficiency Threshold",
    )

    for i, eff in enumerate(avg_efficiency):
        ax2.text(i, eff + 3, f"{eff:.1f}%", ha="center", va="bottom", fontsize=10)

    ax2.set_xlabel("Parallel Configuration", fontsize=12)
    ax2.set_ylabel("Average Parallel Efficiency (%)", fontsize=12)
    ax2.set_title("Efficiency Trend (Large Scale)", fontsize=14, fontweight="bold")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("./image/efficiency_analysis.png", dpi=300, bbox_inches="tight")
    plt.close()


def plot_scalability_analysis():
    """Plot scalability analysis focusing on problem size effects"""
    df = calculate_speedup_efficiency()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Left plot: Strong scalability (fixed problem size, varying processes)
    large_scale = df[df["problem_size"] == "n=131072"]
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
        label="Actual Performance",
    )

    # Ideal scalability
    ideal_times = [avg_times[0] / p for p in processes]
    ax1.plot(
        processes,
        ideal_times,
        "--",
        linewidth=2,
        color="red",
        alpha=0.7,
        label="Ideal Scalability",
    )

    ax1.set_xlabel("Number of Processes", fontsize=12)
    ax1.set_ylabel("Average Execution Time (μs)", fontsize=12)
    ax1.set_title(
        "Strong Scalability Analysis (n=131,072)", fontsize=14, fontweight="bold"
    )
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Right plot: Problem size impact on speedup
    configs = ["1_process", "2_process", "2p_2t", "2p_4t"]
    config_labels = ["1 Process", "2 Processes", "2P+2T", "2P+4T"]

    small_scale = df[df["problem_size"] == "n=4"]
    large_scale = df[df["problem_size"] == "n=131072"]

    small_speedups = []
    large_speedups = []

    for config in configs:
        small_speedup = small_scale[small_scale["config"] == config]["speedup"].mean()
        large_speedup = large_scale[large_scale["config"] == config]["speedup"].mean()
        small_speedups.append(small_speedup)
        large_speedups.append(large_speedup)

    x = np.arange(len(config_labels))
    width = 0.35

    bars1 = ax2.bar(
        x - width / 2,
        small_speedups,
        width,
        label="Small Scale (n=4)",
        alpha=0.8,
        color="lightcoral",
    )
    bars2 = ax2.bar(
        x + width / 2,
        large_speedups,
        width,
        label="Large Scale (n=131,072)",
        alpha=0.8,
        color="steelblue",
    )

    # Add value labels
    for bar, speedup in zip(bars1, small_speedups):
        if speedup < 20:  # Only show reasonable values
            ax2.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.5,
                f"{speedup:.1f}x",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    for bar, speedup in zip(bars2, large_speedups):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.05,
            f"{speedup:.2f}x",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    ax2.set_xlabel("Parallel Configuration", fontsize=12)
    ax2.set_ylabel("Speedup", fontsize=12)
    ax2.set_title("Problem Size Impact on Speedup", fontsize=14, fontweight="bold")
    ax2.set_xticks(x)
    ax2.set_xticklabels(config_labels)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("./image/scalability_analysis.png", dpi=300, bbox_inches="tight")
    plt.close()


def create_performance_summary_table():
    """Create performance summary table"""
    df = calculate_speedup_efficiency()

    # Generate summary table
    summary = (
        df.groupby(["problem_size", "config"])
        .agg({"time": "mean", "speedup": "mean", "efficiency": "mean"})
        .round(2)
    )

    print("Performance Test Results Summary:")
    print("=" * 80)
    print(summary.to_string())
    print("=" * 80)

    # Save as CSV
    summary.to_csv("./image/performance_summary.csv")

    return summary


def main():
    """Main function: generate all charts"""
    print("Generating performance analysis charts...")

    # Generate charts
    plot_speedup_comparison()
    print("✓ Speedup comparison chart saved: ../image/speedup_comparison.png")

    plot_efficiency_analysis()
    print("✓ Efficiency analysis chart saved: ../image/efficiency_analysis.png")

    plot_scalability_analysis()
    print("✓ Scalability analysis chart saved: ../image/scalability_analysis.png")

    # Generate summary table
    create_performance_summary_table()
    print("✓ Performance summary table saved: ../image/performance_summary.csv")

    print("\nAll charts and data tables have been generated successfully!")


if __name__ == "__main__":
    main()
