"""Inputs for reading: files packaged plain, gzip or zip, given as each input kind."""

import codecs
import gzip
import io
import tempfile
import zipfile
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path

import pandas as pd

from snps.utils import gzip_file, zip_file

# packagings of a file at a path: "gzip-fname" as snps.utils.gzip_file writes it,
# with a file name in the gzip header; "zip-stored" as snps.utils.zip_file writes
# it; "-no-ext" names the file without .gz or .zip
PATH_PACKAGINGS = (
    "plain",
    "gzip-fname",
    "gzip-fname-no-ext",
    "zip-stored",
    "zip-deflated-no-ext",
)
# each input kind with each of its packagings; bytes and file objects have no
# name, and between them read a stored and a deflated zip entry
INPUTS = (
    *(("str", packaging) for packaging in PATH_PACKAGINGS),
    *(("path", packaging) for packaging in PATH_PACKAGINGS),
    *(("bytes", packaging) for packaging in ("plain", "gzip", "zip-stored")),
    *(("file", packaging) for packaging in ("plain", "gzip", "zip-deflated")),
)
INPUT_IDS = [f"{kind}-{packaging}" for kind, packaging in INPUTS]


def zip_bytes(
    entries: Mapping[str, bytes], compression: int = zipfile.ZIP_DEFLATED
) -> bytes:
    """Return a zip archive holding ``entries``, in order."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression) as archive:
        for name, data in entries.items():
            archive.writestr(name, data)
    return buffer.getvalue()


def package(plain: Path, packaging: str, directory: Path) -> Path:
    """Return the path of a plain file packaged as ``packaging`` in ``directory``."""
    if packaging == "plain":
        return plain
    directory.mkdir(exist_ok=True)
    extension = ".gz" if packaging.startswith("gzip") else ".zip"
    if packaging.endswith("-no-ext"):
        extension = ""
    path = directory / f"{plain.name}{extension}"
    if packaging.startswith("gzip-fname"):
        gzip_file(str(plain), str(path))
    elif packaging.startswith("gzip"):
        path.write_bytes(gzip.compress(plain.read_bytes()))
    elif packaging.startswith("zip-stored"):
        zip_file(str(plain), str(path), plain.name)
    else:
        path.write_bytes(zip_bytes({plain.name: plain.read_bytes()}))
    return path


@contextmanager
def given_as(path: Path, kind: str) -> Iterator[object]:
    """Give a file as a str path, a pathlib.Path, bytes or an open binary file."""
    if kind == "str":
        yield str(path)
    elif kind == "path":
        yield path
    elif kind == "bytes":
        yield path.read_bytes()
    else:
        with path.open("rb") as f:
            yield f


# file objects that read text: io.TextIOBase streams, and others whose read
# returns str
TEXT_READERS = ("string-io", "text-file", "stream-reader", "spooled-text", "str-read")


class StrReader:
    """File object, not an io.TextIOBase, whose read returns str."""

    closed = False

    def __init__(self, text: str) -> None:
        self._text = text

    def read(self, size=-1):
        return self._text


@contextmanager
def text_reader(kind: str, text: str) -> Iterator[object]:
    """Give text as one of the kinds of file objects that read text."""
    if kind == "string-io":
        yield io.StringIO(text)
    elif kind == "stream-reader":
        yield codecs.getreader("utf-8")(io.BytesIO(text.encode()))
    elif kind == "str-read":
        yield StrReader(text)
    else:
        opener = tempfile.TemporaryFile
        if kind == "spooled-text":
            opener = tempfile.SpooledTemporaryFile
        with opener(mode="w+", encoding="utf-8") as f:
            f.write(text)
            f.seek(0)
            yield f


class NonSeekable(io.RawIOBase):
    """Binary stream of some bytes that cannot seek or tell."""

    def __init__(self, data: bytes) -> None:
        super().__init__()
        self._data = io.BytesIO(data)

    def readable(self):
        return True

    def seekable(self):
        return False

    def seek(self, offset, whence=io.SEEK_SET):
        raise io.UnsupportedOperation("seek")

    def tell(self):
        raise io.UnsupportedOperation("tell")

    def readinto(self, buffer):
        return self._data.readinto(buffer)


def assert_same_data(data, expected) -> None:
    """Assert that two ``SNPData`` hold the same SNPs, phase, source and build."""
    pd.testing.assert_frame_equal(data.frame, expected.frame, check_exact=True)
    pd.testing.assert_series_equal(data.phased, expected.phased, check_exact=True)
    assert data.source == expected.source
    assert data.build == expected.build
    assert data.build_detected == expected.build_detected
