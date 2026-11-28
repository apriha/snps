"""Private, rate-limited JSON transport for resource lookups."""

from __future__ import annotations

import errno
import json
import math
import re
import time
import urllib.error
import urllib.request
from collections import deque
from datetime import timezone
from email.utils import parsedate_to_datetime
from http.client import HTTPException
from threading import Lock
from typing import Never

from snps.exceptions import ResourceError

_ABBREVIATED_WEEKDAY = r"(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun)"
_FULL_WEEKDAY = r"(?:Monday|Tuesday|Wednesday|Thursday|Friday|Saturday|Sunday)"
_MONTH = r"(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)"
_TIME_OF_DAY = r"\d{2}:\d{2}:\d{2}"
_HTTP_DATE = re.compile(
    # IMF-fixdate
    rf"{_ABBREVIATED_WEEKDAY}, \d{{2}} {_MONTH} \d{{4}} {_TIME_OF_DAY} GMT"
    # RFC 850
    rf"|{_FULL_WEEKDAY}, \d{{2}}-{_MONTH}-\d{{2}} {_TIME_OF_DAY} GMT"
    # asctime
    rf"|{_ABBREVIATED_WEEKDAY} {_MONTH} (?:\d{{2}}| \d) {_TIME_OF_DAY} \d{{4}}",
    re.ASCII,
)
_DELAY_SECONDS = re.compile(r"[0-9]+(?:\.[0-9]+)?")
_MAX_RETRIES = 3
_RETRY_FALLBACK = 1.0


def _retry_delay(header: str | None) -> float:
    """Convert a Retry-After field to a delay in seconds.

    References
    ----------
    1. Fielding, Nottingham, and Reschke (2022), HTTP Semantics,
       RFC 9110 §§5.6.7 and 10.2.3.
       https://www.rfc-editor.org/rfc/rfc9110.html#section-10.2.3
    2. Yates (2014), Ensembl REST Rate Limits (fractional seconds).
       https://github.com/Ensembl/ensembl-rest/wiki/Rate-Limits
    """
    if header is None:
        return _RETRY_FALLBACK

    value = header.strip(" \t")
    try:
        if _DELAY_SECONDS.fullmatch(value):
            delay = float(value)
        elif _HTTP_DATE.fullmatch(value):
            date = parsedate_to_datetime(value)
            if date.tzinfo is None:
                date = date.replace(tzinfo=timezone.utc)
            delay = date.timestamp() - time.time()
        else:
            return _RETRY_FALLBACK
    except (ValueError, OverflowError, OSError):
        return _RETRY_FALLBACK

    return max(0.0, delay) if math.isfinite(delay) else _RETRY_FALLBACK


def _sleep_before_retry(delay: float) -> None:
    if delay <= 0:
        return
    try:
        time.sleep(delay)
    except (OverflowError, OSError) as error:
        if delay == _RETRY_FALLBACK or (
            isinstance(error, OSError)
            and error.errno not in (errno.EINVAL, errno.EOVERFLOW)
        ):
            raise
        time.sleep(_RETRY_FALLBACK)


def _reject_json_constant(literal: str) -> Never:
    raise ValueError(f"Invalid JSON constant {literal}")


class RestClient:
    """Fetch JSON objects from one server at a per-instance request rate.

    Parameters
    ----------
    server
        Base URL to which each endpoint is appended.
    requests_per_second
        Maximum number of request starts in a rolling second.
    """

    def __init__(self, *, server: str, requests_per_second: int) -> None:
        if requests_per_second <= 0:
            raise ValueError("requests_per_second must be positive")
        self._server = server
        self._requests_per_second = requests_per_second
        self._starts: deque[float] = deque()
        self._lock = Lock()

    def __getstate__(self) -> tuple[str, int, tuple[float, ...]]:
        with self._lock:
            return self._server, self._requests_per_second, tuple(self._starts)

    def __setstate__(self, state: tuple[str, int, tuple[float, ...]]) -> None:
        self._server, self._requests_per_second, starts = state
        now = time.monotonic()
        self._starts = deque(start for start in starts if start <= now)
        self._lock = Lock()

    def _admit_attempt(self) -> None:
        while True:
            with self._lock:
                now = time.monotonic()
                while self._starts and self._starts[0] + 1.0 <= now:
                    self._starts.popleft()
                if len(self._starts) < self._requests_per_second:
                    self._starts.append(now)
                    return
                delay = self._starts[0] + 1.0 - now
            time.sleep(delay)

    def get_json(self, endpoint: str) -> dict[str, object]:
        """Return the entire JSON object from a GET request.

        Raises
        ------
        ResourceError
            The request failed or its response was unusable.
        """
        url = self._server + endpoint
        for attempt in range(_MAX_RETRIES + 1):
            try:
                request = urllib.request.Request(
                    url, headers={"Content-Type": "application/json"}, method="GET"
                )
                self._admit_attempt()
                with urllib.request.urlopen(request, timeout=30) as response:
                    body = response.read()
            except urllib.error.HTTPError as error:
                status = error.code
                message = f"GET {url} failed with HTTP status {status}"
                try:
                    error.close()
                except (OSError, HTTPException, ValueError):
                    raise ResourceError(message) from error
                if status != 429 or attempt == _MAX_RETRIES:
                    raise ResourceError(message) from error
                retry_after = (
                    error.headers.get("Retry-After") if error.headers else None
                )
                delay = _retry_delay(retry_after)
                try:
                    _sleep_before_retry(delay)
                except (OSError, OverflowError) as sleep_error:
                    raise ResourceError(message) from sleep_error
            except (OSError, HTTPException, ValueError, OverflowError) as error:
                raise ResourceError(f"GET {url} failed") from error
            else:
                try:
                    result: object = json.loads(
                        body.decode("utf-8"), parse_constant=_reject_json_constant
                    )
                    if not isinstance(result, dict):
                        raise ValueError("Expected a JSON object")
                except (ValueError, RecursionError) as error:
                    raise ResourceError(
                        f"GET {url} returned an unusable response"
                    ) from error
                return result
        raise AssertionError("Request attempt loop ended without a result")
