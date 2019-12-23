import logging
import os
import random

from atomicwrites import atomic_write
import pandas as pd

from snps import SNPs
from snps.resources import Resources
from snps.utils import Parallelizer, save_df_as_csv, create_dir

OUTPUT_DIR = "output"

# create output directory for this example
create_dir(OUTPUT_DIR)

# assume script is being run from examples dir
r = Resources(resources_dir="../../resources")

# setup logger to output to file in output directory
logging.basicConfig(
    filename="{}".format(os.path.join(OUTPUT_DIR, "xy-chrom-snp-ratios.txt")),
    format="%(asctime)s - %(message)s",
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger()


def get_xy_chrom_snp_ratios(task):
    file = task["file"]

    try:
        s = SNPs(r.load_opensnp_datadump_file(file), assign_par_snps=False)
    except Exception as err:
        logger.error("{}: {}".format(file, err))
        return None

    try:
        if s.snp_count != 0:
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
                s.source,
                s.build,
                s.build_detected,
                s.determine_sex(),
                x_snps,
                heterozygous_x_snps,
                y_snps,
                y_snps_not_null,
                s.snp_count,
                s.chromosomes_summary,
                file,
            ]
        else:
            logger.info("{}: {}".format(file, "no SNPs processed"))

    except Exception as err:
        logger.error("{}: {}".format(file, err))
        return None


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
            "source",
            "build",
            "build_detected",
            "sex",
            "x_snps",
            "heterozygous_x_snps",
            "y_snps",
            "y_snps_not_null",
            "snp_count",
            "chromosomes_summary",
            "file",
        ],
    )

    # derive the columns we want to analyze
    df["heterozygous_x_snps_ratio"] = df.heterozygous_x_snps / df.x_snps
    df["y_snps_not_null_ratio"] = df.y_snps_not_null / df.y_snps

    # create the histograms
    hist = df.hist(
        column=["heterozygous_x_snps_ratio", "y_snps_not_null_ratio"],
        grid=False,
        bins=15,
        figsize=(8, 6),
        edgecolor="black",
    )
    hist[0, 0].set_ylabel("Frequency")
    hist[0, 1].set_ylabel("Frequency")
    fig = hist[0, 0].get_figure()
    fig.tight_layout()

    # save output
    with atomic_write(
        "{}".format(os.path.join(OUTPUT_DIR, "xy-chrom-snp-ratios.png")),
        mode="wb",
        overwrite=True,
    ) as f:
        fig.savefig(f)
    save_df_as_csv(df, OUTPUT_DIR, "xy-chrom-snp-ratios.csv")

    logger.info("stop")
