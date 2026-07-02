"""Pytest configuration: offline resources by default + README doctests via Sybil.

An autouse fixture injects a fixture-backed resource provider so the entire suite
(including the README examples) runs with zero network and zero mocks. The Sybil
configuration enables testing Python code blocks in README.md as part of pytest:
``PythonCodeBlockParser`` evaluates fenced ``python`` blocks, while ``SkipParser``
allows selectively skipping examples with Markdown comments.
"""

import pytest
from sybil import Sybil
from sybil.parsers.markdown import PythonCodeBlockParser, SkipParser

from snps.resources import set_default_provider
from tests.support import FakeResources


@pytest.fixture(autouse=True)
def _offline_resources():
    """Use a fixture-backed resource provider for every test (no network, no mocks)."""
    set_default_provider(FakeResources())
    try:
        yield
    finally:
        set_default_provider(None)


pytest_collect_file = Sybil(
    parsers=[
        PythonCodeBlockParser(),
        SkipParser(),
    ],
    patterns=["README.md"],
    fixtures=["_offline_resources"],
).pytest()
