""" Parse openSNP datadump files.

Attempt to parse each file in the openSNP datadump. For files where SNPs were loaded,
save summary statistics to a dataframe and output as a CSV. For files where no SNPs were
loaded, save a message for each file indicating the issue and optionally extract these
files from the datadump for debugging.
"""

import logging
import os
import random

from atomicwrites import atomic_write
import pandas as pd

from snps import SNPs
from snps.resources import Resources
from snps.utils import Parallelizer, save_df_as_csv, create_dir, clean_str

OUTPUT_DIR = "output"
EXTRACT_FILES = True

# create output directory for this example
create_dir(OUTPUT_DIR)

# assume script is being run from examples dir
r = Resources(resources_dir="../../resources")

# setup logger to output to file in output directory
logging.basicConfig(
    filename="{}".format(os.path.join(OUTPUT_DIR, "parse-opensnp-files.txt")),
    format="%(asctime)s: %(message)s",
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger()


def load_file(task):
    file = task["file"]

    try:
        s = SNPs(r.load_opensnp_datadump_file(file), assign_par_snps=False)
    except Exception as err:
        return {"msg": str(err).strip()[:100], "file": file}

    if s.snp_count != 0:
        d = s.get_summary()
        d.update({"file": file})
        return d
    else:
        return {"msg": "no SNPs processed", "file": file}


def main():
    logger.info("start")

    # get filenames from openSNP data dump
    filenames = r.get_opensnp_datadump_filenames()

    filenames = [
        filename
        for filename in filenames
        if "readme" not in filename and "phenotype" not in filename
    ]

    # draw a sample from the observations
    random.seed(1)
    SAMPLE_SIZE = len(filenames)
    # SAMPLE_SIZE = 10
    samples = random.sample(range(len(filenames)), SAMPLE_SIZE)

    # setup tasks for parallelizing / execution on multiple cores
    p = Parallelizer(parallelize=True)
    tasks = [{"file": filenames[i]} for i in samples]

    # run tasks; results is a list of dicts
    results = p(load_file, tasks)

    # get results from `load_file` where `snp_count` was non-zero
    rows = [item for item in results if "msg" not in item]

    df = pd.DataFrame(
        rows,
        columns=[
            "file",
            "source",
            "build",
            "build_detected",
            "chromosomes",
            "snp_count",
        ],
    )

    save_df_as_csv(df, OUTPUT_DIR, "parse-opensnp-files.csv")

    # log parsing statistics
    file_count = len(filenames)
    logger.info("{} files in the openSNP datadump".format(file_count))
    logger.info("{:.2%} of openSNP datadump files parsed".format(len(df) / file_count))
    logger.info(
        "build detected in {:.2%} of files parsed".format(
            len(df.loc[df.build_detected]) / len(df)
        )
    )

    # extract files from the datadump where `load_file` returned a message
    if EXTRACT_FILES:
        # group files with same message (e.g., {"some message": ["file1", "file2"], ...})
        d = {}
        for result in results:
            if "msg" in result:
                if result["msg"] in d:
                    d[result["msg"]].append(result["file"])
                else:
                    d[result["msg"]] = [result["file"]]

        # add messages / file filters as necessary...
        d["build not detected"] = list(df.loc[~df.build_detected].file.values)

        # extract files that have messages for debugging
        for msg, files in d.items():
            if len(files) == 0:
                continue

            # create a directory for each message (prefix indicates number of files)
            path = os.path.join(
                OUTPUT_DIR, "{:04}_{}".format(len(files), clean_str(msg))
            )
            create_dir(path)
            # save each file with message into created directory
            for filename in files:
                with atomic_write(os.path.join(path, filename), mode="wb") as f:
                    f.write(r.load_opensnp_datadump_file(filename))

    logger.info("stop")


if __name__ == "__main__":
    main()
