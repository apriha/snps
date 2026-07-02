"""Fixture-backed ``ResourceProvider`` for offline tests (no network, no mocks)."""

import atexit
import os
import shutil
import tempfile

from snps.resources import Resources
from snps.utils import create_dir, gzip_file
from tests.support.data import PAR_EFFECTS, standard_assembly_mappings

_TESTS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_RESOURCES_DIR = os.path.join(_TESTS_DIR, "resources")
_INPUT_DIR = os.path.join(_TESTS_DIR, "input")

# committed raw fixtures served in place of downloads, keyed by cache filename
_FIXTURES = {
    "gsa_rsid_map.txt.gz": "gsa_rsid_map.txt",
    "gsa_chrpos_map.txt.gz": "gsa_chrpos_map.txt",
    "dbsnp_151_37_reverse.txt.gz": "dbsnp_151_37_reverse.txt",
    "chip_clusters.tsv.gz": "chip_clusters.tsv",
    "low_quality_snps.tsv.gz": "low_quality_snps.tsv",
}

# shared cache dir so committed fixtures are gzipped once and reused across instances
_CACHE_DIR = tempfile.mkdtemp(prefix="snps-fake-resources-")
atexit.register(shutil.rmtree, _CACHE_DIR, ignore_errors=True)

_PAR_LOOKUPS = {e["refsnp_id"]: e for e in PAR_EFFECTS}


class FakeResources(Resources):
    """Offline ``ResourceProvider`` backed by committed fixtures.

    Inherits all parsing / transforms / aggregators from :class:`~snps.resources.Resources`
    and overrides only the byte source (``_fetch``) plus the two REST-built methods
    (``get_assembly_mapping_data``, ``get_par_lookup``). Optional keyword arguments inject
    pre-built data for tests that need specific resources.

    Parameters
    ----------
    chip_clusters : pandas.DataFrame, optional
        value returned by ``get_chip_clusters`` (defaults to the committed fixture)
    low_quality_snps : pandas.DataFrame, optional
        value returned by ``get_low_quality_snps`` (defaults to the committed fixture)
    assembly_mapping_data : dict or callable, optional
        value returned by ``get_assembly_mapping_data``; a callable is invoked with
        ``(source_assembly, target_assembly)`` (defaults to the standard test mappings)
    resources_dir : str, optional
        cache directory used when gzipping committed fixtures (defaults to a shared tmp dir)
    """

    def __init__(
        self,
        *,
        chip_clusters=None,
        low_quality_snps=None,
        assembly_mapping_data=None,
        resources_dir=None,
    ):
        super().__init__(resources_dir=resources_dir or _CACHE_DIR)
        self._chip_clusters = chip_clusters
        self._low_quality_snps = low_quality_snps
        self._assembly_mapping_data = assembly_mapping_data

    def _fetch(self, url, filename):
        base = os.path.basename(filename)
        if base in _FIXTURES:
            src = os.path.join(_RESOURCES_DIR, _FIXTURES[base])
        elif filename.replace(os.sep, "/").startswith("fasta/") and base.endswith(
            ".fa.gz"
        ):
            src = os.path.join(_INPUT_DIR, "generic.fa")
        else:
            return ""

        destination = os.path.join(self._resources_dir, filename)
        create_dir(os.path.dirname(destination))
        if not os.path.exists(destination):
            gzip_file(src, destination)
        return destination

    def get_assembly_mapping_data(self, source_assembly, target_assembly):
        if self._assembly_mapping_data is not None:
            data = self._assembly_mapping_data
            if callable(data):
                return data(source_assembly, target_assembly)
            return data
        return standard_assembly_mappings().get((source_assembly, target_assembly), {})

    def get_par_lookup(self, rsid):
        snapshot = _PAR_LOOKUPS.get(rsid.split("rs")[1])
        if snapshot is None:
            return None
        if "merged_snapshot_data" in snapshot:
            merged = "rs" + snapshot["merged_snapshot_data"]["merged_into"][0]
            return self.get_par_lookup(merged)
        if "nosnppos_snapshot_data" in snapshot:
            return None
        return snapshot
