import logging
import os
import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from atomicwrites import atomic_write
from matplotlib import patches
from snps import SNPs
from snps.resources import Resources
from snps.utils import Parallelizer, create_dir, save_df_as_csv

OUTPUT_DIR = "output"

# create output directory for this example
create_dir(OUTPUT_DIR)

# assume script is being run from examples dir
r = Resources(resources_dir="../../resources")

# setup logger to output to file in output directory
logging.basicConfig(
    filename=f"{os.path.join(OUTPUT_DIR, 'xy-chrom-snp-ratios.txt')}",
    format="%(asctime)s#%(message)s",
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger()


def get_xy_chrom_snp_ratios(task):
    file = task["file"]

    try:
        logger.info(f"loading {file}")
        s = SNPs(r.load_opensnp_datadump_file(file), assign_par_snps=False)
    except Exception as err:
        logger.error(f"{file}#{err}")
        return None

    try:
        if s.count != 0:
            # get X chromosome statistics
            x_snps = len(s.snps.loc[(s.snps["chrom"] == "X")])
            heterozygous_x_snps = len(
                s.snps.loc[
                    (s.snps["chrom"] == "X")
                    & (s.snps["genotype"].notnull())
                    & (s.snps["genotype"].str.len() == 2)
                    & (s.snps["genotype"].str[0] != s.snps["genotype"].str[1])
                ]
            )

            # get Y chromosome statistics
            y_snps = len(s.snps.loc[(s.snps["chrom"] == "Y")])
            y_snps_not_null = len(
                s.snps.loc[(s.snps["chrom"] == "Y") & (s.snps["genotype"].notnull())]
            )

            return [
                file,
                s.source,
                s.build,
                s.build_detected,
                x_snps,
                heterozygous_x_snps,
                y_snps,
                y_snps_not_null,
                s.count,
            ]
        else:
            logger.info(f"{file}#{'no SNPs processed'}")

    except Exception as err:
        logger.error(f"{file}#{err}")
        return None


def create_analysis_plot(
    df_ratios, heterozygous_x_snps_threshold=0.03, y_snps_not_null_threshold=0.3
):
    # https://matplotlib.org/gallery/lines_bars_and_markers/scatter_hist.html#sphx-glr-gallery-lines-bars-and-markers-scatter-hist-py

    x = df_ratios["heterozygous_x_snps_ratio"]
    y = df_ratios["y_snps_not_null_ratio"]
    alpha = 0.1

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # start with a rectangular Figure
    fig = plt.figure(figsize=(8, 8))
    fig.suptitle(f"Analysis of openSNP datadump XY chrom SNP ratios; N = {len(df)}")

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction="in", top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction="in", labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction="in", labelleft=False)

    # now determine nice limits by hand:
    bins = 40
    ax_spacing = 0.025
    binwidth_x = x.max() / bins
    binwidth_y = y.max() / bins
    lim_x = np.ceil(x.max() / binwidth_x) * binwidth_x
    lim_y = np.ceil(y.max() / binwidth_y) * binwidth_y
    ax_scatter.set_xlim((-lim_x * ax_spacing, lim_x + lim_x * ax_spacing))
    ax_scatter.set_ylim((-lim_y * ax_spacing, lim_y + lim_y * ax_spacing))
    ax_scatter.grid(True)
    ax_scatter.set_axisbelow(True)
    ax_scatter.set_xlabel("Heterozygous X SNPs / X SNPs")
    ax_scatter.set_ylabel("Y SNPs not null / Y SNPs")

    heterozygous_x_snps_threshold_line = ax_scatter.axvline(
        x=heterozygous_x_snps_threshold,
        c="blue",
        label=f"Het. X threshold={heterozygous_x_snps_threshold}",
    )
    y_snps_not_null_threshold_line = ax_scatter.axhline(
        y=y_snps_not_null_threshold,
        c="red",
        label=f"Y not null threshold={y_snps_not_null_threshold}",
    )

    # fill genotype areas
    ax_scatter.fill_between(
        [heterozygous_x_snps_threshold, lim_x + lim_x * ax_spacing],
        y_snps_not_null_threshold,
        facecolor="blue",
        alpha=alpha,
    )
    ax_scatter.fill_between(
        [0, heterozygous_x_snps_threshold],
        y_snps_not_null_threshold,
        lim_y + lim_y * ax_spacing,
        facecolor="red",
        alpha=alpha,
    )

    # add the points
    ax_scatter.scatter(x, y, s=6)

    # add the histograms
    bins_x = np.arange(0, lim_x + binwidth_x, binwidth_x)
    bins_y = np.arange(0, lim_y + binwidth_y, binwidth_y)
    ax_histx.hist(x, bins=bins_x, edgecolor="black")
    ax_histy.hist(y, bins=bins_y, orientation="horizontal", edgecolor="black")
    ax_histx.axvline(x=heterozygous_x_snps_threshold, c="blue")
    ax_histy.axhline(y=y_snps_not_null_threshold, c="red")
    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())

    # add legend
    handles = [
        heterozygous_x_snps_threshold_line,
        y_snps_not_null_threshold_line,
        patches.Patch(color="blue", alpha=alpha, label="Female genotypes"),
        patches.Patch(color="red", alpha=alpha, label="Male genotypes"),
    ]
    fig.legend(
        handles=handles, loc="upper right", bbox_to_anchor=(0.99, 0.95), fontsize=8
    )

    # annotate count of values in each quadrant of the scatterplot
    x_offset = lim_x * 0.01
    y_offset = lim_y * 0.01
    ax_scatter.annotate(
        f"n={len(df_ratios.loc[(df_ratios.heterozygous_x_snps_ratio < heterozygous_x_snps_threshold) & (df_ratios.y_snps_not_null_ratio < y_snps_not_null_threshold)])}",
        (
            heterozygous_x_snps_threshold - x_offset,
            y_snps_not_null_threshold - y_offset,
        ),
        ha="right",
        va="top",
    )
    ax_scatter.annotate(
        f"n={ len(df_ratios.loc[(df_ratios.heterozygous_x_snps_ratio < heterozygous_x_snps_threshold)& (df_ratios.y_snps_not_null_ratio >= y_snps_not_null_threshold)])}",
        (
            heterozygous_x_snps_threshold - x_offset,
            y_snps_not_null_threshold + y_offset,
        ),
        ha="right",
        va="bottom",
    )
    ax_scatter.annotate(
        f"n={len(df_ratios.loc[ (df_ratios.heterozygous_x_snps_ratio >= heterozygous_x_snps_threshold) & (df_ratios.y_snps_not_null_ratio >= y_snps_not_null_threshold)])}",
        (
            heterozygous_x_snps_threshold + x_offset,
            y_snps_not_null_threshold + y_offset,
        ),
        ha="left",
        va="bottom",
    )
    ax_scatter.annotate(
        f"n={len(df_ratios.loc[(df_ratios.heterozygous_x_snps_ratio >= heterozygous_x_snps_threshold) & (df_ratios.y_snps_not_null_ratio < y_snps_not_null_threshold)])}",
        (
            heterozygous_x_snps_threshold + x_offset,
            y_snps_not_null_threshold - y_offset,
        ),
        ha="left",
        va="top",
    )

    return plt


if __name__ == "__main__":
    logger.info("start")

    # get filenames from openSNP data dump
    filenames = r.get_opensnp_datadump_filenames()

    # draw a sample from the observations
    random.seed(1)
    SAMPLE_SIZE = len(filenames)
    # SAMPLE_SIZE = 10
    samples = random.sample(range(len(filenames)), SAMPLE_SIZE)

    # setup tasks for parallelizing / execution on multiple cores
    p = Parallelizer(parallelize=True)
    tasks = [{"file": filenames[i]} for i in samples]

    # results are a list of lists
    rows = p(get_xy_chrom_snp_ratios, tasks)

    # remove None results
    rows = [row for row in rows if row]

    df = pd.DataFrame(
        rows,
        columns=[
            "file",
            "source",
            "build",
            "build_detected",
            "x_snps",
            "heterozygous_x_snps",
            "y_snps",
            "y_snps_not_null",
            "count",
        ],
    )

    # derive the columns we want to analyze
    df["heterozygous_x_snps_ratio"] = df.heterozygous_x_snps / df.x_snps
    df["y_snps_not_null_ratio"] = df.y_snps_not_null / df.y_snps

    df.drop(df.loc[df["heterozygous_x_snps_ratio"].isna()].index, inplace=True)
    df.drop(df.loc[df["y_snps_not_null_ratio"].isna()].index, inplace=True)

    plt = create_analysis_plot(
        df[["heterozygous_x_snps_ratio", "y_snps_not_null_ratio"]]
    )

    # save output
    with atomic_write(
        f"{os.path.join(OUTPUT_DIR, 'xy-chrom-snp-ratios.png')}",
        mode="wb",
        overwrite=True,
    ) as f:
        plt.savefig(f)

    save_df_as_csv(df, OUTPUT_DIR, "xy-chrom-snp-ratios.csv")

    logger.info("stop")
