#!/usr/bin/env python3

import json
import sys
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt


def load_results(path: str):
    with open(path) as f:
        return json.load(f)


def build_axes_values(results, x_axis: str, y_axis: str):
    xs = set()
    ys = set()
    for row in results["configuration"]:
        xs.add(row[x_axis])
        ys.add(row[y_axis])
    xs = sorted(list(xs))
    ys = sorted(list(ys))
    return xs, ys


def build_grid(results, x_axis: str, y_axis: str, key: str, xs, ys):
    Z = np.full((len(ys), len(xs)), np.nan)
    for row in results["configuration"]:
        try:
            val = row["Result"][key]
        except KeyError:
            continue
        Z[ys.index(row[y_axis]), xs.index(row[x_axis])] = val
    return Z


def metric_title(metric: str) -> str:
    mapping = {
        "APY": "APY",
        "APY_boost": "APY - boost",
        "slippage": "Slippage",
        "volume": "Volume",
    }
    return mapping.get(metric, metric.replace("_", " ").title())


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Plot selected metrics as an Nx3 grid (max 3 per row)")
    parser.add_argument(
        "metrics",
        nargs="*",
        help="Metric keys under results['Result'] (e.g., APY APY_boost slippage volume)",
    )
    parser.add_argument(
        "-f",
        "--file",
        default="results.json",
        help="Path to results.json (default: results.json)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)

    metrics = args.metrics or [ "xcp_profit", "donation_coin0_total", "n_rebalances", "trade_volume", "n_trades", "cex_follow_time_frac"]
    results = load_results(args.file)

    x_axis = "A"
    y_axis = "mid_fee"

    xs, ys = build_axes_values(results, x_axis, y_axis)

    n = len(metrics)
    cols = min(3, max(1, n))
    rows = int(math.ceil(n / 3.0))
    # Scale figure size with rows/cols
    fig_width = 4.2 * cols
    fig_height = 3.6 * rows
    fig, axes_grid = plt.subplots(rows, cols, figsize=(fig_width, fig_height), squeeze=False)
    axes = axes_grid.flatten()
    cmap = plt.get_cmap("jet")

    for i, metric in enumerate(metrics):
        ax = axes[i]
        Z = build_grid(results, x_axis, y_axis, metric, xs, ys)

        ax.set_xscale("log")
        ax.set_yscale("log")
        im = ax.pcolormesh(xs, ys, Z, cmap=cmap)
        im.set_edgecolor("face")
        cbar = fig.colorbar(im, ax=ax)
        ax.set_xlabel(x_axis)
        ax.set_ylabel(y_axis)
        title = metric_title(metric)
        cbar.set_label(title, rotation=270, labelpad=15)
        ax.set_title(title)

    # Hide any unused axes (if metrics not multiple of cols)
    for j in range(n, rows * cols):
        axes[j].axis('off')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
