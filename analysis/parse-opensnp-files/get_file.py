"""Get a file from the openSNP datadump for debugging."""

import os

from atomicwrites import atomic_write

from snps.resources import Resources
from snps.utils import create_dir

OUTPUT_DIR = "output"
FILE = "user662_file340_yearofbirth_unknown_sex_unknown.23andme.txt"

if __name__ == "__main__":
    # create output directory for this example
    create_dir(OUTPUT_DIR)

    # assume script is being run from examples dir
    r = Resources(resources_dir="../../resources")

    with atomic_write(os.path.join(OUTPUT_DIR, FILE), mode="wb") as f:
        f.write(r.load_opensnp_datadump_file(FILE))
