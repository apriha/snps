import errno
import io
import json
import pickle
import urllib.error
from collections import deque
from copy import deepcopy
from email.message import Message
from email.utils import formatdate

import pytest

from snps import _rest
from snps.exceptions import ResourceError, SNPsError
from snps.resources import Resources

ENSEMBL = "https://rest.ensembl.org"
NCBI = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp"
ENDPOINT = "/synthetic/rs990001"
URL = ENSEMBL + ENDPOINT
WALL_TIME = 1_700_000_000


class FakeTransport:
    def __init__(self):
        self.pending = deque()
        self.calls = []
        self.starts = []
        self.sleeps = []
        self.now = 0.0
        self.wall_time = WALL_TIME
        self.failed_bodies = []

    @property
    def urls(self):
        return [request.full_url for request, _ in self.calls]

    def urlopen(self, request, timeout):
        self.calls.append((request, timeout))
        self.starts.append(self.now)
        assert self.pending, f"Unexpected request to {request.full_url}"
        response = self.pending.popleft()
        if isinstance(response, urllib.error.HTTPError):
            self.failed_bodies.append(response.fp)
        if isinstance(response, Exception):
            raise response
        return response

    def monotonic(self):
        return self.now

    def time(self):
        return self.wall_time + self.now

    def sleep(self, seconds):
        assert all(body.closed for body in self.failed_bodies)
        self.sleeps.append(seconds)
        self.now += seconds

    def queue_bytes(self, data):
        body = io.BytesIO(data)
        self.pending.append(body)
        return body

    def queue_json(self, data):
        return self.queue_bytes(json.dumps(data, ensure_ascii=False).encode("utf-8"))

    def queue_http_error(self, url, status, retry_after=None):
        headers = Message()
        if retry_after is not None:
            headers["Retry-After"] = retry_after
        error = urllib.error.HTTPError(
            url, status, "Synthetic HTTP error", headers, io.BytesIO(b"error")
        )
        self.pending.append(error)
        return error


@pytest.fixture
def transport(monkeypatch):
    fake = FakeTransport()
    monkeypatch.setattr(_rest.urllib.request, "urlopen", fake.urlopen)
    monkeypatch.setattr(_rest, "time", fake)
    return fake


@pytest.fixture
def rest_client(transport):
    return _rest.RestClient(server=ENSEMBL, requests_per_second=15)


def test_resource_error_directly_subclasses_snps_error():
    assert ResourceError.__bases__ == (SNPsError,)


def test_get_sends_json_header_and_decodes_utf8_object(transport):
    server = "https://transport.invalid"
    endpoint = "/synthetic/a%2Fb/rs990001"
    result = {"id": "rs990001", "label": "café", "nested": {"count": 2}}
    body = transport.queue_json(result)

    client = _rest.RestClient(server=server, requests_per_second=15)
    assert client.get_json(endpoint) == result

    [(request, timeout)] = transport.calls
    assert request.get_method() == "GET"
    assert timeout == 30
    assert request.full_url == server + endpoint
    assert [
        value
        for name, value in request.header_items()
        if name.lower() == "content-type"
    ] == ["application/json"]
    assert body.closed


def test_http_429_recovers_after_fractional_delay_and_reuses_url(
    transport, rest_client
):
    error = transport.queue_http_error(URL, 429, "0.125")
    body = transport.queue_json({})

    assert rest_client.get_json(ENDPOINT) == {}

    assert transport.urls == [URL, URL]
    assert [timeout for _, timeout in transport.calls] == [30, 30]
    assert transport.sleeps == [0.125]
    assert error.fp.closed
    assert body.closed


def test_http_429_close_failure_aborts_retry(transport, rest_client, monkeypatch):
    error = transport.queue_http_error(URL, 429, "0.125")

    def fail_close():
        raise OSError("Synthetic close failure")

    with error.fp as body, monkeypatch.context() as patch:
        patch.setattr(body, "close", fail_close)
        with pytest.raises(ResourceError) as raised:
            rest_client.get_json(ENDPOINT)

    assert URL in str(raised.value)
    assert "429" in str(raised.value)
    assert raised.value.__cause__ is error
    assert transport.urls == [URL]
    assert transport.sleeps == []


def test_http_429_sleep_range_error_falls_back_to_one_second(
    transport, rest_client, monkeypatch
):
    delay = 9223372036
    transport.queue_http_error(URL, 429, str(delay))
    transport.queue_json({})
    sleep = transport.sleep
    rejected_delays = []

    def sleep_with_limit(seconds):
        if seconds == delay:
            rejected_delays.append(seconds)
            raise OSError(errno.EINVAL, "Synthetic unrepresentable timeout")
        sleep(seconds)

    monkeypatch.setattr(transport, "sleep", sleep_with_limit)

    assert rest_client.get_json(ENDPOINT) == {}

    assert rejected_delays == [delay]
    assert transport.sleeps == [1]
    assert transport.urls == [URL, URL]
    assert transport.starts == [0, 1]


@pytest.mark.parametrize(
    "retry_after",
    [
        None,
        "invalid",
        "-0.25",
        "NaN",
        "Infinity",
        "1_0",
        "１２",
        "+45",
        "1e2",
        formatdate(WALL_TIME + 3, usegmt=True) + ", junk",
    ],
)
def test_http_429_exhausts_three_retries_with_one_second_fallback(
    transport, rest_client, retry_after
):
    errors = [transport.queue_http_error(URL, 429, retry_after) for _ in range(4)]

    with pytest.raises(ResourceError) as raised:
        rest_client.get_json(ENDPOINT)

    assert URL in str(raised.value)
    assert "429" in str(raised.value)
    assert raised.value.__cause__ is errors[-1]
    assert transport.urls == [URL] * 4
    assert transport.sleeps == [1, 1, 1]
    assert all(error.fp.closed for error in errors)


@pytest.mark.parametrize(
    ("wall_time", "retry_after", "delay"),
    [
        (WALL_TIME, formatdate(WALL_TIME + 3, usegmt=True), 3),
        (WALL_TIME, formatdate(WALL_TIME - 3, usegmt=True), 0),
        # Each legacy date represents 1994-11-06 08:49:37 UTC, three seconds later.
        (784_111_774, "Sunday, 06-Nov-94 08:49:37 GMT", 3),
        (784_111_774, "Sun Nov  6 08:49:37 1994", 3),
        (784_111_774, "Sun Nov 06 08:49:37 1994", 3),
    ],
)
def test_http_429_retry_after_date_uses_wall_time(
    transport, rest_client, wall_time, retry_after, delay
):
    transport.wall_time = wall_time
    transport.queue_http_error(URL, 429, retry_after)
    transport.queue_json({"id": "rs990001"})

    assert rest_client.get_json(ENDPOINT) == {"id": "rs990001"}

    assert transport.starts == [0, delay]


def test_http_503_is_terminal_and_its_start_counts_toward_limit(transport, capsys):
    error = transport.queue_http_error(URL, 503)
    client = _rest.RestClient(server=ENSEMBL, requests_per_second=1)

    with pytest.raises(ResourceError) as raised:
        client.get_json(ENDPOINT)

    assert URL in str(raised.value)
    assert "503" in str(raised.value)
    assert raised.value.__cause__ is error
    assert transport.urls == [URL]
    assert transport.sleeps == []
    assert error.fp.closed
    assert capsys.readouterr().err == ""

    transport.queue_json({})
    assert client.get_json(ENDPOINT) == {}
    assert transport.starts == [0, 1]


@pytest.mark.parametrize("error_type", [urllib.error.URLError, TimeoutError])
def test_transport_error_names_url_and_preserves_cause(
    transport, rest_client, error_type
):
    error = error_type("Synthetic transport failure")
    transport.pending.append(error)

    with pytest.raises(ResourceError) as raised:
        rest_client.get_json(ENDPOINT)

    assert URL in str(raised.value)
    assert raised.value.__cause__ is error
    assert transport.urls == [URL]


@pytest.mark.parametrize(
    ("data", "cause_type"),
    [
        (b"", None),
        (b"{", json.JSONDecodeError),
        (b"\xff", UnicodeDecodeError),
        (b"[]", None),
        (b"null", None),
    ],
)
def test_invalid_response_raises_resource_error(
    transport, rest_client, data, cause_type
):
    body = transport.queue_bytes(data)

    with pytest.raises(ResourceError) as raised:
        rest_client.get_json(ENDPOINT)

    assert URL in str(raised.value)
    if cause_type is not None:
        assert isinstance(raised.value.__cause__, cause_type)
    assert transport.urls == [URL]
    assert body.closed


def test_rolling_limit_counts_retries_and_is_per_client(transport):
    client = _rest.RestClient(server=ENSEMBL, requests_per_second=2)
    transport.queue_json({})
    transport.queue_http_error(URL, 429, "0.125")
    for _ in range(3):
        transport.queue_json({})

    assert client.get_json(ENDPOINT) == {}
    transport.now += 0.25
    transport.wall_time += 3600  # Wall time can jump without aging the rate window.
    assert client.get_json(ENDPOINT) == {}
    assert client.get_json(ENDPOINT) == {}
    other_client = _rest.RestClient(server=ENSEMBL, requests_per_second=2)
    assert other_client.get_json(ENDPOINT) == {}

    # Starts at 0 and 0.25 release capacity at 1 and 1.25 respectively.
    assert transport.starts == [0, 0.25, 1, 1.25, 1.25]


def test_pickled_rest_client_allows_request_after_monotonic_reset(transport):
    transport.now = 86_410
    client = _rest.RestClient(server=ENSEMBL, requests_per_second=1)
    transport.queue_json({"value": 1})
    transport.queue_json({"value": 2})

    assert client.get_json(ENDPOINT) == {"value": 1}
    state = pickle.dumps(client)
    transport.now = 10
    restored = pickle.loads(state)

    assert restored.get_json(ENDPOINT) == {"value": 2}
    assert transport.starts == [86_410, 10]
    assert transport.sleeps == []


def test_refsnp_follows_first_merge_id_and_paces_later_lookups(transport, tmp_path):
    transport.queue_json(
        {"merged_snapshot_data": {"merged_into": ["990003", "990002"]}}
    )
    transport.queue_json(
        {"merged_snapshot_data": {"merged_into": ["990004", "990005"]}}
    )
    primary = {
        "refsnp_id": "990004",
        "primary_snapshot_data": {"placements_with_allele": []},
        "opaque": {"example": ["preserved"]},
    }
    transport.queue_json(primary)
    subsequent = {**primary, "refsnp_id": "990007"}
    transport.queue_json({"nosnppos_snapshot_data": {"example": "unplaced"}})
    transport.queue_json(subsequent)
    resources = Resources(resources_dir=tmp_path)

    assert resources.get_par_lookup("rs990001") == primary
    assert resources.get_par_lookup("rs990006") is None
    assert resources.get_par_lookup("rs990007") == subsequent

    assert transport.urls == [
        f"{NCBI}/{rsid}" for rsid in ("990001", "990003", "990004", "990006", "990007")
    ]
    # All five starts share the NCBI limit of one request per second.
    assert transport.starts == [0, 1, 2, 3, 4]


@pytest.mark.parametrize("after_merge", [False, True])
def test_refsnp_transport_failure_raises_resource_error(
    transport, tmp_path, after_merge
):
    failed_id = "990003" if after_merge else "990001"
    if after_merge:
        transport.queue_json({"merged_snapshot_data": {"merged_into": ["990003"]}})
    failed_url = f"{NCBI}/{failed_id}"
    error = transport.queue_http_error(failed_url, 503)
    resources = Resources(resources_dir=tmp_path)

    with pytest.raises(ResourceError) as raised:
        resources.get_par_lookup("rs990001")

    assert failed_url in str(raised.value)
    expected_urls = [f"{NCBI}/990001", failed_url] if after_merge else [failed_url]
    assert transport.urls == expected_urls
    assert error.fp.closed


@pytest.mark.parametrize(
    ("merged_into", "message"),
    [([], "Invalid RefSNP merge target"), (["990001"], "RefSNP merge cycle")],
)
def test_refsnp_invalid_merge_stops_after_one_lookup(
    transport, tmp_path, merged_into, message
):
    transport.queue_json({"merged_snapshot_data": {"merged_into": merged_into}})
    resources = Resources(resources_dir=tmp_path)

    with pytest.raises(ResourceError, match=message):
        resources.get_par_lookup("rs990001")

    assert transport.urls == [f"{NCBI}/990001"]


def test_resources_can_be_deepcopied_and_pickled_after_lookup(transport, tmp_path):
    resources = Resources(resources_dir=tmp_path)
    unplaced = {"nosnppos_snapshot_data": {"example": "unplaced"}}
    transport.queue_json(unplaced)
    assert resources.get_par_lookup("rs990001") is None

    copied = deepcopy(resources)
    restored = pickle.loads(pickle.dumps(resources))
    primary = {
        "refsnp_id": "990002",
        "primary_snapshot_data": {"placements_with_allele": []},
    }
    transport.queue_json(primary)
    transport.queue_json(unplaced)

    assert copied.get_par_lookup("rs990002") == primary
    assert restored.get_par_lookup("rs990003") is None
    assert transport.urls == [
        f"{NCBI}/{rsid}" for rsid in ("990001", "990002", "990003")
    ]


def test_assembly_mapping_rolls_back_then_caches_all_chromosomes(transport, tmp_path):
    chromosomes = [str(number) for number in range(1, 23)] + ["X", "Y", "MT"]
    urls = [
        f"{ENSEMBL}/map/human/NCBI36/{chromosome}/GRCh37" for chromosome in chromosomes
    ]
    payloads = [{"mappings": [], "opaque": chromosome} for chromosome in chromosomes]
    archive = tmp_path / "NCBI36_GRCh37.tar.gz"
    resources = Resources(resources_dir=tmp_path)
    initial_files = {path for path in tmp_path.rglob("*") if path.is_file()}
    transport.queue_json(payloads[0])
    error = transport.queue_http_error(urls[1], 503)

    with pytest.raises(ResourceError) as raised:
        resources.get_assembly_mapping_data("NCBI36", "GRCh37")

    assert "503" in str(raised.value)
    assert urls[1] in str(raised.value)
    assert [url.removesuffix("?") for url in transport.urls] == urls[:2]
    assert error.fp.closed
    assert not archive.exists()
    assert {path for path in tmp_path.rglob("*") if path.is_file()} == initial_files

    for payload in payloads:
        transport.queue_json(payload)
    result = resources.get_assembly_mapping_data("NCBI36", "GRCh37")

    assert [url.removesuffix("?") for url in transport.urls] == urls[:2] + urls
    assert result == dict(zip(chromosomes, payloads))
    assert list(result) == chromosomes
    assert archive.is_file()
    # Two starts from the failed run leave thirteen slots in the first second.
    assert transport.starts == [0] * 15 + [1] * 12

    cached_resources = Resources(resources_dir=tmp_path)
    assert cached_resources.get_assembly_mapping_data("NCBI36", "GRCh37") == result
    assert len(transport.calls) == 27
