# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Ensembl REST client.

Notes
-----
Modified from https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client.

References
----------
1. Yates et. al. (doi:10.1093/bioinformatics/btu613),
   `<http://europepmc.org/search/?query=DOI:10.1093/bioinformatics/btu613>`_
2. Zerbino et. al. (doi.org/10.1093/nar/gkx1098), https://doi.org/10.1093/nar/gkx1098

"""

import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


class EnsemblRestClient:
    def __init__(self, server="https://rest.ensembl.org", reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if "Content-Type" not in hdrs:
            hdrs["Content-Type"] = "application/json"

        if params:
            endpoint += "?" + urllib.parse.urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = urllib.request.Request(self.server + endpoint, headers=hdrs)

            with urllib.request.urlopen(request) as response:
                content = response.read().decode("utf-8")

            if content:
                data = json.loads(content)
            self.req_count += 1

        except urllib.error.HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if "Retry-After" in e.headers:
                    retry = e.headers["Retry-After"]
                    time.sleep(float(retry))
                    return self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write(
                    f"Request failed for {endpoint}: Status code: {e.code} Reason: {e.reason}\n"
                )

        return data
