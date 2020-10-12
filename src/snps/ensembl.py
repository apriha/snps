""" Ensembl REST client.

Notes
-----
Modified from https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client.

References
----------
1. Yates et. al. (doi:10.1093/bioinformatics/btu613),
   `<http://europepmc.org/search/?query=DOI:10.1093/bioinformatics/btu613>`_
2. Zerbino et. al. (doi.org/10.1093/nar/gkx1098), https://doi.org/10.1093/nar/gkx1098

"""

"""
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and

"""

"""
BSD 3-Clause License

Copyright (c) 2019, Andrew Riha
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write(
                    f"Request failed for {endpoint}: Status code: {e.code} Reason: {e.reason}\n"
                )

        return data
