#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# 设置字体
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def load_data(filename):
    """Load performance test data"""
    df = pd.read_csv(filename)
    return df


def create_performance_comparison(df):
    """Create performance comparison charts"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle("NTT Unified Framework Performance Test Results", fontsize=16)

    # Get all algorithms and parallel strategies
    algorithms = df["Algorithm"].unique()
    parallel_strategies = df["Parallel_Strategy"].unique()
    sizes = sorted(df["Size"].unique())

    # Create subplot for each algorithm
    for i, algorithm in enumerate(algorithms):
        row = i // 2
        col = i % 2
        ax = axes[row, col]

        algorithm_data = df[df["Algorithm"] == algorithm]

        for strategy in parallel_strategies:
            strategy_data = algorithm_data[
                algorithm_data["Parallel_Strategy"] == strategy
            ]
            if not strategy_data.empty:
                ax.plot(
                    strategy_data["Size"],
                    strategy_data["Avg_Time_ms"],
                    marker="o",
                    label=strategy,
                    linewidth=2,
                    markersize=6,
                )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Polynomial Length")
        ax.set_ylabel("Average Time (ms)")
        ax.set_title(f"{algorithm} Algorithm Performance Comparison")
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("performance_comparison.png", dpi=300, bbox_inches="tight")
    plt.show()


def create_speedup_analysis(df):
    """Create speedup analysis charts"""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("Parallel Speedup Analysis", fontsize=16)

    algorithms = df["Algorithm"].unique()

    for i, algorithm in enumerate(algorithms):
        ax = axes[i]
        algorithm_data = df[df["Algorithm"] == algorithm]

        # Get serial version as baseline
        serial_data = algorithm_data[algorithm_data["Parallel_Strategy"] == "Serial"]

        if serial_data.empty:
            continue

        # Calculate speedup
        for strategy in ["OpenMP", "GPU"]:
            strategy_data = algorithm_data[
                algorithm_data["Parallel_Strategy"] == strategy
            ]
            if not strategy_data.empty:
                # Merge data to calculate speedup
                merged = pd.merge(
                    serial_data,
                    strategy_data,
                    on="Size",
                    suffixes=("_serial", "_parallel"),
                )
                speedup = merged["Avg_Time_ms_serial"] / merged["Avg_Time_ms_parallel"]

                ax.plot(
                    merged["Size"],
                    speedup,
                    marker="o",
                    label=f"{strategy} vs Serial",
                    linewidth=2,
                    markersize=6,
                )

        ax.set_xscale("log")
        ax.set_xlabel("Polynomial Length")
        ax.set_ylabel("Speedup")
        ax.set_title(f"{algorithm} Algorithm Speedup")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axhline(y=1, color="red", linestyle="--", alpha=0.5, label="No Speedup")

    plt.tight_layout()
    plt.savefig("speedup_analysis.png", dpi=300, bbox_inches="tight")
    plt.show()


def create_heatmap(df):
    """Create performance heatmap"""
    # Select data with maximum size for heatmap analysis
    max_size = df["Size"].max()
    heatmap_data = df[df["Size"] == max_size]

    # Create pivot table
    pivot_data = heatmap_data.pivot(
        index="Algorithm", columns="Parallel_Strategy", values="Avg_Time_ms"
    )

    plt.figure(figsize=(10, 6))
    sns.heatmap(
        pivot_data,
        annot=True,
        fmt=".2f",
        cmap="YlOrRd",
        cbar_kws={"label": "Average Time (ms)"},
    )
    plt.title(f"Performance Heatmap for Polynomial Length {max_size}")
    plt.tight_layout()
    plt.savefig("performance_heatmap.png", dpi=300, bbox_inches="tight")
    plt.show()


def create_strategy_grouped_comparison(df):
    """Create performance comparison charts grouped by parallel strategy"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle("NTT Performance Analysis by Parallel Strategy", fontsize=16)

    # Get all algorithms and parallel strategies
    algorithms = df["Algorithm"].unique()
    parallel_strategies = df["Parallel_Strategy"].unique()
    sizes = sorted(df["Size"].unique())

    # Create subplot for each parallel strategy
    for i, strategy in enumerate(parallel_strategies):
        row = i // 2
        col = i % 2
        ax = axes[row, col]

        strategy_data = df[df["Parallel_Strategy"] == strategy]

        for algorithm in algorithms:
            algorithm_data = strategy_data[strategy_data["Algorithm"] == algorithm]
            if not algorithm_data.empty:
                ax.plot(
                    algorithm_data["Size"],
                    algorithm_data["Avg_Time_ms"],
                    marker="o",
                    label=algorithm,
                    linewidth=2,
                    markersize=6,
                )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Polynomial Length")
        ax.set_ylabel("Average Time (ms)")
        ax.set_title(f"{strategy} Strategy Performance Comparison")
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("strategy_grouped_comparison.png", dpi=300, bbox_inches="tight")
    plt.show()


def print_summary_statistics(df):
    """Print summary statistics"""
    print("=== NTT Unified Framework Performance Test Summary ===")
    print()

    # Group statistics by algorithm and parallel strategy
    summary = (
        df.groupby(["Algorithm", "Parallel_Strategy"])
        .agg({"Avg_Time_ms": ["mean", "min", "max"], "Size": ["min", "max"]})
        .round(3)
    )

    print("Performance Statistics Summary:")
    print(summary)
    print()

    # Calculate best performance combination
    best_performance = df.loc[df["Avg_Time_ms"].idxmin()]
    print(f"Best Performance Combination:")
    print(f"  Algorithm: {best_performance['Algorithm']}")
    print(f"  Parallel Strategy: {best_performance['Parallel_Strategy']}")
    print(f"  Size: {best_performance['Size']}")
    print(f"  Average Time: {best_performance['Avg_Time_ms']:.3f} ms")
    print()

    # Calculate speedup ratios
    for algorithm in df["Algorithm"].unique():
        algorithm_data = df[df["Algorithm"] == algorithm]
        serial_times = algorithm_data[algorithm_data["Parallel_Strategy"] == "Serial"][
            "Avg_Time_ms"
        ]

        for strategy in ["OpenMP", "GPU"]:
            strategy_times = algorithm_data[
                algorithm_data["Parallel_Strategy"] == strategy
            ]["Avg_Time_ms"]
            if not strategy_times.empty and not serial_times.empty:
                avg_speedup = serial_times.mean() / strategy_times.mean()
                print(f"{algorithm} + {strategy} Average Speedup: {avg_speedup:.2f}x")


def main():
    """Main function"""
    try:
        # Load data
        df = load_data("performance_results.csv")
        print("Data loaded successfully!")
        print(f"Total tests: {len(df)}")
        print()

        # Print summary statistics
        print_summary_statistics(df)

        # Create algorithm-grouped charts (original)
        print("Generating algorithm-grouped performance comparison charts...")
        create_performance_comparison(df)

        # Create strategy-grouped charts (new)
        print("Generating strategy-grouped performance comparison charts...")
        create_strategy_grouped_comparison(df)

        print("Generating speedup analysis charts...")
        create_speedup_analysis(df)

        print("Generating performance heatmap...")
        create_heatmap(df)

        print("Analysis completed! Charts saved as PNG files.")

    except Exception as e:
        print(f"Error: {e}")
        print(
            "Please ensure performance_results.csv file exists and has correct format."
        )


if __name__ == "__main__":
    main()
