"""Client for NCBI Variation Services RefSNP snapshots."""

from __future__ import annotations

import logging
from collections.abc import Mapping

from snps._rest import RestClient
from snps.exceptions import ResourceError

logger = logging.getLogger(__name__)
_SERVER = "https://api.ncbi.nlm.nih.gov"
_REFSNP_PATH = "/variation/v0/refsnp/"


class VariationClient:
    """Look up RefSNP snapshots through one rate-limited NCBI client."""

    def __init__(self) -> None:
        self._rest_client = RestClient(server=_SERVER, requests_per_second=1)

    def get_par_lookup(self, rsid: str) -> Mapping[str, object] | None:
        """Resolve a RefSNP snapshot for a PAR SNP."""
        seen: set[str] = set()
        while True:
            number = rsid.split("rs")[1]
            endpoint = _REFSNP_PATH + number
            if number in seen:
                raise ResourceError(f"RefSNP merge cycle at {_SERVER}{endpoint}")
            seen.add(number)

            snapshot = self._rest_client.get_json(endpoint)
            if "merged_snapshot_data" in snapshot:
                metadata = snapshot["merged_snapshot_data"]
                targets: object = (
                    metadata.get("merged_into") if isinstance(metadata, dict) else None
                )
                if (
                    not isinstance(targets, list)
                    or not targets
                    or not isinstance(targets[0], str)
                    or not targets[0].isascii()
                    or not targets[0].isdecimal()
                ):
                    raise ResourceError(
                        f"Invalid RefSNP merge target from {_SERVER}{endpoint}"
                    )
                merged_rsid = "rs" + targets[0]
                logger.info("SNP id %s has been merged into id %s", rsid, merged_rsid)
                rsid = merged_rsid
            elif "nosnppos_snapshot_data" in snapshot:
                logger.warning("Unable to look up SNP id %s", rsid)
                return None
            else:
                return snapshot
