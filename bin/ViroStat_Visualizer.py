#!/usr/bin/env python3
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import LogLocator, NullFormatter

def plot_viral_raw_data(input_file: str) -> None:
    if not os.path.exists(input_file):
        print(f"Error: 找不到文件 '{input_file}'", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(input_file, sep='\t')

    required_cols = {"Groupname", "Count"}
    if not required_cols.issubset(df.columns):
        print("Error: 文件必须包含 'Groupname' 和 'Count' 列", file=sys.stderr)
        sys.exit(1)

    if df.empty:
        print("Error: 输入文件为空，无法作图", file=sys.stderr)
        sys.exit(1)

    sample_counts = df["Groupname"].value_counts()
    order = sample_counts.index.tolist()
    labels = [f"{name}\n(n={sample_counts[name]})" for name in order]

    sns.set_theme(style="whitegrid")
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial"]
    plt.rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(figsize=(14, 8))

    sns.boxplot(
        data=df,
        x="Groupname",
        y="Count",
        order=order,
        ax=ax,
        hue="Groupname",
        palette="Set3",
        legend=False,
        fliersize=0,
        width=0.6,
        linewidth=1.5
    )

    sns.stripplot(
        data=df,
        x="Groupname",
        y="Count",
        order=order,
        ax=ax,
        color=".3",
        size=3,
        alpha=0.5,
        jitter=True
    )

    ax.set_yscale("log")
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs="auto", numticks=12))
    ax.yaxis.set_minor_formatter(NullFormatter())

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")

    ax.set_xlabel("")
    ax.set_ylabel("Read Counts (Log$_{10}$ Scale)", fontsize=12, fontweight="bold")

    output_png = os.path.join(os.path.dirname(os.path.abspath(input_file)), "Viral_Distribution_Final_Clean.png")
    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    print(f"Success! 图像已保存为: {output_png}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python ViroStat_Visualizer.py Groupname_Sample_Count_Detail.xls", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    plot_viral_raw_data(input_file)

if __name__ == "__main__":
    main()
    