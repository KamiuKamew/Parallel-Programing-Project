#!/usr/bin/env python3
"""
Performance Analysis and Visualization Script for MPI NTT Experiment
Handles the vast difference in speedup values between different problem sizes
by using separate charts and appropriate scaling techniques.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set matplotlib style
plt.style.use("default")
plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["axes.facecolor"] = "white"
plt.rcParams["font.size"] = 10


def calculate_speedup_efficiency():
    """Calculate speedup and efficiency from performance data"""
    # Performance data from manual testing
    performance_data = [
        # Small scale (n=4)
        {"problem_size": "n=4", "config": "serial", "time": 2.0, "modulus": "p1"},
        {"problem_size": "n=4", "config": "1_process", "time": 0.0119, "modulus": "p1"},
        {"problem_size": "n=4", "config": "2_process", "time": 0.0132, "modulus": "p1"},
        {"problem_size": "n=4", "config": "2p_2t", "time": 0.0169, "modulus": "p1"},
        {"problem_size": "n=4", "config": "2p_4t", "time": 0.1387, "modulus": "p1"},
        # Large scale (n=131072)
        {
            "problem_size": "n=131072",
            "config": "serial",
            "time": 473.0,
            "modulus": "p1",
        },
        {
            "problem_size": "n=131072",
            "config": "1_process",
            "time": 489.0,
            "modulus": "p1",
        },
        {
            "problem_size": "n=131072",
            "config": "2_process",
            "time": 277.0,
            "modulus": "p1",
        },
        {"problem_size": "n=131072", "config": "2p_2t", "time": 274.0, "modulus": "p1"},
        {"problem_size": "n=131072", "config": "2p_4t", "time": 311.0, "modulus": "p1"},
    ]

    df = pd.DataFrame(performance_data)

    # Calculate speedup and efficiency for each problem size
    results = []

    for size in df["problem_size"].unique():
        size_data = df[df["problem_size"] == size]
        serial_time = size_data[size_data["config"] == "serial"]["time"].iloc[0]

        for _, row in size_data.iterrows():
            if row["config"] != "serial":
                speedup = serial_time / row["time"]

                # Calculate theoretical process count for efficiency
                if row["config"] == "1_process":
                    processes = 1
                elif row["config"] == "2_process":
                    processes = 2
                elif row["config"] == "2p_2t":
                    processes = 4  # 2 processes × 2 threads
                elif row["config"] == "2p_4t":
                    processes = 8  # 2 processes × 4 threads

                efficiency = (speedup / processes) * 100

                results.append(
                    {
                        "problem_size": row["problem_size"],
                        "config": row["config"],
                        "time": row["time"],
                        "speedup": speedup,
                        "efficiency": efficiency,
                    }
                )

    return pd.DataFrame(results)


def plot_separate_scale_analysis():
    """Plot speedup analysis with separate charts for different problem scales"""
    df = calculate_speedup_efficiency()

    # Create a 2x2 subplot layout
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    configs = ["1_process", "2_process", "2p_2t", "2p_4t"]
    config_labels = [
        "1 Process\n(MPI)",
        "2 Processes\n(MPI)",
        "2P+2T\n(Hybrid)",
        "2P+4T\n(Hybrid)",
    ]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    # Top Left: Small Scale Speedup (Linear Scale)
    small_scale = df[df["problem_size"] == "n=4"]
    small_speedups = []

    for config in configs:
        speedup = small_scale[small_scale["config"] == config]["speedup"].iloc[0]
        small_speedups.append(speedup)

    bars1 = ax1.bar(
        config_labels,
        small_speedups,
        color=colors,
        alpha=0.8,
        edgecolor="black",
        linewidth=0.5,
    )

    # Add value labels
    for bar, speedup in zip(bars1, small_speedups):
        ax1.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(small_speedups) * 0.02,
            f"{speedup:.0f}x",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
        )

    ax1.set_ylabel("Speedup Factor", fontsize=12)
    ax1.set_title(
        "Small Scale Performance (n=4)\nLinear Scale", fontsize=14, fontweight="bold"
    )
    ax1.grid(True, alpha=0.3, axis="y")
    ax1.set_ylim(0, max(small_speedups) * 1.15)

    # Top Right: Large Scale Speedup (Linear Scale)
    large_scale = df[df["problem_size"] == "n=131072"]
    large_speedups = []

    for config in configs:
        speedup = large_scale[large_scale["config"] == config]["speedup"].iloc[0]
        large_speedups.append(speedup)

    bars2 = ax2.bar(
        config_labels,
        large_speedups,
        color=colors,
        alpha=0.8,
        edgecolor="black",
        linewidth=0.5,
    )

    # Add ideal speedup line
    ideal_processes = [1, 2, 4, 8]
    ax2.plot(
        range(len(config_labels)),
        ideal_processes,
        "r--",
        linewidth=2,
        alpha=0.7,
        label="Theoretical Maximum",
        marker="o",
        markersize=4,
    )

    # Add value labels
    for bar, speedup in zip(bars2, large_speedups):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.05,
            f"{speedup:.2f}x",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
        )

    ax2.set_ylabel("Speedup Factor", fontsize=12)
    ax2.set_title(
        "Large Scale Performance (n=131,072)\nLinear Scale",
        fontsize=14,
        fontweight="bold",
    )
    ax2.grid(True, alpha=0.3, axis="y")
    ax2.legend(loc="upper right")
    ax2.set_ylim(0, max(max(large_speedups), max(ideal_processes)) * 1.15)

    # Bottom Left: Combined Comparison (Log Scale)
    x = np.arange(len(config_labels))
    width = 0.35

    bars3 = ax3.bar(
        x - width / 2,
        small_speedups,
        width,
        label="Small Scale (n=4)",
        alpha=0.8,
        color="lightcoral",
        edgecolor="black",
        linewidth=0.5,
    )
    bars4 = ax3.bar(
        x + width / 2,
        large_speedups,
        width,
        label="Large Scale (n=131,072)",
        alpha=0.8,
        color="steelblue",
        edgecolor="black",
        linewidth=0.5,
    )

    ax3.set_yscale("log")
    ax3.set_xlabel("Parallel Configuration", fontsize=12)
    ax3.set_ylabel("Speedup Factor (Log Scale)", fontsize=12)
    ax3.set_title(
        "Combined Scale Comparison\nLogarithmic Scale", fontsize=14, fontweight="bold"
    )
    ax3.set_xticks(x)
    ax3.set_xticklabels(config_labels)
    ax3.legend()
    ax3.grid(True, alpha=0.3, which="both")

    # Add horizontal reference lines on log scale
    ax3.axhline(y=1, color="gray", linestyle=":", alpha=0.5, label="No Speedup")
    ax3.axhline(y=10, color="orange", linestyle=":", alpha=0.5, label="10x Speedup")
    ax3.axhline(y=100, color="red", linestyle=":", alpha=0.5, label="100x Speedup")

    # Bottom Right: Efficiency Comparison
    small_efficiency = []
    large_efficiency = []

    for config in configs:
        small_eff = small_scale[small_scale["config"] == config]["efficiency"].iloc[0]
        large_eff = large_scale[large_scale["config"] == config]["efficiency"].iloc[0]
        small_efficiency.append(small_eff)
        large_efficiency.append(large_eff)

    bars5 = ax4.bar(
        x - width / 2,
        small_efficiency,
        width,
        label="Small Scale (n=4)",
        alpha=0.8,
        color="lightcoral",
        edgecolor="black",
        linewidth=0.5,
    )
    bars6 = ax4.bar(
        x + width / 2,
        large_efficiency,
        width,
        label="Large Scale (n=131,072)",
        alpha=0.8,
        color="steelblue",
        edgecolor="black",
        linewidth=0.5,
    )

    # Add efficiency reference lines
    ax4.axhline(y=100, color="r", linestyle="--", alpha=0.7, label="Perfect Efficiency")
    ax4.axhline(
        y=80, color="orange", linestyle="--", alpha=0.7, label="Good Efficiency"
    )
    ax4.axhline(
        y=50, color="gray", linestyle="--", alpha=0.7, label="Acceptable Efficiency"
    )

    ax4.set_xlabel("Parallel Configuration", fontsize=12)
    ax4.set_ylabel("Parallel Efficiency (%)", fontsize=12)
    ax4.set_title(
        "Efficiency Comparison\nAcross Problem Scales", fontsize=14, fontweight="bold"
    )
    ax4.set_xticks(x)
    ax4.set_xticklabels(config_labels)
    ax4.legend(loc="upper right")
    ax4.grid(True, alpha=0.3, axis="y")

    plt.tight_layout(pad=2.0)
    plt.savefig("./image/separate_scale_analysis.png", dpi=300, bbox_inches="tight")
    plt.close()

    print("Created separate scale analysis chart: ./image/separate_scale_analysis.png")


def plot_detailed_performance_breakdown():
    """Create detailed performance breakdown focusing on meaningful comparisons"""
    df = calculate_speedup_efficiency()

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 10))

    configs = ["1_process", "2_process", "2p_2t", "2p_4t"]
    config_labels = ["1 Process", "2 Processes", "2P+2T", "2P+4T"]

    # 1. Small Scale Only (n=4) - Detailed Analysis
    small_scale = df[df["problem_size"] == "n=4"]
    small_times = [2.0]  # Serial baseline
    small_speedups = [1.0]  # Serial baseline

    for config in configs:
        time = small_scale[small_scale["config"] == config]["time"].iloc[0]
        speedup = small_scale[small_scale["config"] == config]["speedup"].iloc[0]
        small_times.append(time)
        small_speedups.append(speedup)

    labels_with_serial = ["Serial"] + config_labels
    bars1 = ax1.bar(
        labels_with_serial,
        small_speedups,
        color=["gray"] + ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"],
        alpha=0.8,
        edgecolor="black",
        linewidth=0.5,
    )

    for i, (bar, speedup) in enumerate(zip(bars1, small_speedups)):
        if i == 0:  # Serial
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 2,
                f"{speedup:.1f}x",
                ha="center",
                va="bottom",
                fontsize=10,
            )
        else:
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 5,
                f"{speedup:.0f}x",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )

    ax1.set_ylabel("Speedup Factor", fontsize=12)
    ax1.set_title(
        "Small Scale (n=4): All Configurations", fontsize=14, fontweight="bold"
    )
    ax1.grid(True, alpha=0.3, axis="y")
    ax1.tick_params(axis="x", rotation=45)

    # 2. Large Scale Only (n=131072) - Detailed Analysis
    large_scale = df[df["problem_size"] == "n=131072"]
    large_times = [473.0]  # Serial baseline
    large_speedups = [1.0]  # Serial baseline

    for config in configs:
        time = large_scale[large_scale["config"] == config]["time"].iloc[0]
        speedup = large_scale[large_scale["config"] == config]["speedup"].iloc[0]
        large_times.append(time)
        large_speedups.append(speedup)

    bars2 = ax2.bar(
        labels_with_serial,
        large_speedups,
        color=["gray"] + ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"],
        alpha=0.8,
        edgecolor="black",
        linewidth=0.5,
    )

    # Add theoretical maximum
    theoretical = [1, 1, 2, 4, 8]
    ax2.plot(
        range(len(labels_with_serial)),
        theoretical,
        "r--",
        linewidth=2,
        alpha=0.7,
        label="Theoretical Max",
        marker="o",
        markersize=4,
    )

    for bar, speedup in zip(bars2, large_speedups):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.05,
            f"{speedup:.2f}x",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    ax2.set_ylabel("Speedup Factor", fontsize=12)
    ax2.set_title(
        "Large Scale (n=131,072): All Configurations", fontsize=14, fontweight="bold"
    )
    ax2.grid(True, alpha=0.3, axis="y")
    ax2.legend()
    ax2.tick_params(axis="x", rotation=45)

    # 3. Execution Time Comparison (Log Scale)
    ax3.bar(
        range(len(labels_with_serial)),
        small_times,
        width=0.4,
        label="Small Scale (n=4)",
        alpha=0.8,
        color="lightcoral",
    )
    ax3.bar(
        [x + 0.4 for x in range(len(labels_with_serial))],
        large_times,
        width=0.4,
        label="Large Scale (n=131,072)",
        alpha=0.8,
        color="steelblue",
    )

    ax3.set_yscale("log")
    ax3.set_xlabel("Configuration", fontsize=12)
    ax3.set_ylabel("Execution Time (μs, Log Scale)", fontsize=12)
    ax3.set_title("Execution Time Comparison", fontsize=14, fontweight="bold")
    ax3.set_xticks([x + 0.2 for x in range(len(labels_with_serial))])
    ax3.set_xticklabels(labels_with_serial, rotation=45)
    ax3.legend()
    ax3.grid(True, alpha=0.3, which="both")

    # 4. Problem Size Effect Analysis
    problem_size_effect = []
    for i, config in enumerate(["serial"] + configs):
        if config == "serial":
            ratio = large_times[0] / small_times[0]
        else:
            ratio = large_times[i] / small_times[i]
        problem_size_effect.append(ratio)

    bars4 = ax4.bar(
        labels_with_serial,
        problem_size_effect,
        color=["gray"] + ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"],
        alpha=0.8,
        edgecolor="black",
        linewidth=0.5,
    )

    for bar, ratio in zip(bars4, problem_size_effect):
        ax4.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 500,
            f"{ratio:.0f}x",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    ax4.set_ylabel("Time Ratio (Large/Small)", fontsize=12)
    ax4.set_title("Problem Size Impact Factor", fontsize=14, fontweight="bold")
    ax4.grid(True, alpha=0.3, axis="y")
    ax4.tick_params(axis="x", rotation=45)

    plt.tight_layout(pad=2.0)
    plt.savefig(
        "./image/detailed_performance_breakdown.png", dpi=300, bbox_inches="tight"
    )
    plt.close()

    print(
        "Created detailed performance breakdown: ./image/detailed_performance_breakdown.png"
    )


def create_performance_summary_table():
    """Create a comprehensive performance summary table"""
    df = calculate_speedup_efficiency()

    # Reorganize data for the summary table
    summary_data = []

    for size in ["n=4", "n=131072"]:
        size_data = df[df["problem_size"] == size]
        size_label = "Small Scale (n=4)" if size == "n=4" else "Large Scale (n=131,072)"

        # Add serial baseline
        if size == "n=4":
            serial_time = 2.0
        else:
            serial_time = 473.0

        summary_data.append(
            {
                "Problem Size": size_label,
                "Configuration": "Serial Baseline",
                "Time (μs)": f"{serial_time:.1f}",
                "Speedup": "1.00x",
                "Efficiency": "100.0%",
            }
        )

        # Add parallel configurations
        config_names = {
            "1_process": "1 MPI Process",
            "2_process": "2 MPI Processes",
            "2p_2t": "2 Processes + 2 Threads",
            "2p_4t": "2 Processes + 4 Threads",
        }

        for _, row in size_data.iterrows():
            summary_data.append(
                {
                    "Problem Size": size_label,
                    "Configuration": config_names[row["config"]],
                    "Time (μs)": f"{row['time']:.2f}",
                    "Speedup": f"{row['speedup']:.2f}x",
                    "Efficiency": f"{row['efficiency']:.1f}%",
                }
            )

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv("./performance_summary_separated.csv", index=False)

    print("Created performance summary table: ./performance_summary_separated.csv")
    print("\nPerformance Summary:")
    print(summary_df.to_string(index=False))


def main():
    """Generate all performance analysis charts and tables"""
    print("Generating performance analysis with proper scale separation...")

    # Create output directory if it doesn't exist
    import os

    os.makedirs("./image", exist_ok=True)

    # Generate all charts
    plot_separate_scale_analysis()
    plot_detailed_performance_breakdown()
    create_performance_summary_table()

    print("\n=== Analysis Complete ===")
    print("Generated files:")
    print("- ./image/separate_scale_analysis.png")
    print("- ./image/detailed_performance_breakdown.png")
    print("- ./performance_summary_separated.csv")

    print("\nKey Findings:")
    print(
        "• Small scale problems (n=4) show extremely high speedups due to overhead domination"
    )
    print(
        "• Large scale problems (n=131,072) show realistic parallel performance characteristics"
    )
    print(
        "• 2-process MPI configuration provides optimal balance for large-scale problems"
    )
    print(
        "• Hybrid parallelization shows diminishing returns with increased thread count"
    )


if __name__ == "__main__":
    main()
