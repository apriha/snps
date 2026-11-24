import gzip
import io
import lzma
import struct
import zipfile
import zlib
from functools import partial

import pytest

import snps.io.readers._source
from snps import SNPs
from snps.exceptions import ParseError, UnsupportedFormatError
from snps.io import Reader
from snps.io.readers import read_file
from snps.io.readers._registry import Header
from snps.io.readers._source import InputSource
from tests.support.cases import CASES, INPUT, assert_matches_case
from tests.support.inputs import NonSeekable, assert_same_data, package, zip_bytes

FIXTURE = INPUT / "23andme.txt"
CASE = next(case for case in CASES if case.id == "23andme")
GFG = "GFG_filtered_unphased_genotypes_23andMe.txt"
UNRECOGNIZED = b"invented notes\nnothing here\n"
PACKAGINGS = ("plain", "gzip", "zip")


@pytest.fixture
def data():
    return FIXTURE.read_bytes()


@pytest.fixture
def expected():
    return read_file(FIXTURE)


def packed(data, packaging):
    """Return data packaged plain, gzip or zip, as bytes."""
    if packaging == "gzip":
        return gzip.compress(data)
    if packaging == "zip":
        return zip_bytes({"a.txt": data})
    return data


def stored_zip(data, *, flags=0, method=None, crc_change=0):
    """Return a zip of one stored entry, with fields of its headers changed."""
    content = bytearray(zip_bytes({"a.txt": data}, zipfile.ZIP_STORED))
    # offsets of the flags, method and CRC fields in the local and central headers
    for signature, offsets in (
        (b"PK\x03\x04", (6, 8, 14)),
        (b"PK\x01\x02", (8, 10, 16)),
    ):
        start = content.index(signature)
        flag_at, method_at, crc_at = (start + offset for offset in offsets)
        (flag_bits,) = struct.unpack_from("<H", content, flag_at)
        struct.pack_into("<H", content, flag_at, flag_bits | flags)
        if method is not None:
            struct.pack_into("<H", content, method_at, method)
        (crc,) = struct.unpack_from("<I", content, crc_at)
        struct.pack_into("<I", content, crc_at, crc ^ crc_change)
    return bytes(content)


def _entry_data_start(content):
    """Return where the data of the first zip entry starts, after its local header."""
    name_length, extra_length = struct.unpack_from("<HH", content, 26)
    return 30 + name_length + extra_length


def _stored_gzip(data):
    """Return a gzip member holding data in one stored (uncompressed) deflate block."""
    # magic, deflate, no flags, no modification time, no extra flags, unknown system
    header = b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff"
    # the final stored block: its type, its length and the complement of its length
    block = b"\x01" + struct.pack("<HH", len(data), len(data) ^ 0xFFFF) + data
    trailer = struct.pack("<II", zlib.crc32(data), len(data))
    member = bytearray(header + block + trailer)
    assert gzip.decompress(member) == data
    return member


def _gzip_corrupt_body(data):
    # a changed data byte breaks the checksum in the trailer
    member = _stored_gzip(data)
    member[-9] ^= 0x01
    return bytes(member)


def _gzip_corrupt_block(data):
    # a length complement that does not match the length is not deflate data
    member = _stored_gzip(data)
    member[13] ^= 0xFF
    return bytes(member)


def _zip_corrupt_bzip2(data):
    # a bzip2 stream starts with "BZh"
    content = bytearray(zip_bytes({"a.txt": data}, zipfile.ZIP_BZIP2))
    start = _entry_data_start(content)
    assert content[start : start + 3] == b"BZh"
    content[start + 2] ^= 0xFF
    return bytes(content)


def _zip_corrupt_lzma(data):
    # a version and the size of the properties come before the properties, whose
    # first byte must be below 225
    content = bytearray(zip_bytes({"a.txt": data}, zipfile.ZIP_LZMA))
    start = _entry_data_start(content)
    assert struct.unpack_from("<H", content, start + 2) == (5,)
    content[start + 4] = 0xFF
    return bytes(content)


def _zip_central_directory_offset(data):
    # the end record closes with the central directory offset and a comment length
    content = bytearray(zip_bytes({"a.txt": data}))
    content[-5] ^= 0x01
    return bytes(content)


def _zip_undecodable_local_name(data):
    # the UTF-8 flag on a name that is not UTF-8, in the local header only
    content = bytearray(zip_bytes({"a.txt": data}))
    (flag_bits,) = struct.unpack_from("<H", content, 6)
    struct.pack_into("<H", content, 6, flag_bits | 0x800)
    content[30:35] = b"\xff\xfe.tx"
    return bytes(content)


# compression, a function making the input from the fixture's data, and the type
# of the error that the corrupt data error is chained from
CORRUPT = {
    # truncated inside its first lines, so that detection finds it too
    "truncated-gzip": ("gzip", lambda data: gzip.compress(data)[:20], EOFError),
    "gzip-corrupt-body": ("gzip", _gzip_corrupt_body, zlib.error),
    "gzip-corrupt-block": ("gzip", _gzip_corrupt_block, zlib.error),
    "gzip-then-bad-member": (
        "gzip",
        lambda data: gzip.compress(data) + b"\x1f\x8binvented junk",
        zlib.error,
    ),
    "zip-signature-only": (
        "zip",
        lambda data: b"PK\x03\x04invented junk",
        zipfile.BadZipFile,
    ),
    "truncated-zip": (
        "zip",
        lambda data: zip_bytes({"a.txt": data})[:-30],
        zipfile.BadZipFile,
    ),
    "zip-bad-crc": (
        "zip",
        partial(stored_zip, crc_change=0xFFFF),
        zipfile.BadZipFile,
    ),
    "zip-corrupt-bzip2": ("zip", _zip_corrupt_bzip2, OSError),
    "zip-corrupt-lzma": ("zip", _zip_corrupt_lzma, lzma.LZMAError),
    "zip-central-directory-offset": ("zip", _zip_central_directory_offset, None),
    "zip-undecodable-local-name": (
        "zip",
        _zip_undecodable_local_name,
        UnicodeDecodeError,
    ),
}


def assert_corrupt_data_error(exc, name, compression, cause):
    assert type(exc) is ParseError
    assert str(exc) == (
        f"Cannot read {name}: the {compression} compressed data is corrupt or "
        "truncated."
    )
    if cause is None:
        assert exc.__cause__ is None
    else:
        assert isinstance(exc.__cause__, cause)


@pytest.mark.parametrize("name", ["x.txt", "x.txt.gz", "x.zip"])
@pytest.mark.parametrize(
    ("packaging", "compression"), [("plain", None), ("gzip", "gzip"), ("zip", "zip")]
)
def test_compression_from_leading_bytes(tmp_path, data, name, packaging, compression):
    path = tmp_path / name
    path.write_bytes(packed(data, packaging))
    with InputSource(path) as source:
        assert source.compression == compression
        with source.open() as stream:
            assert stream.read() == data


@pytest.mark.parametrize(
    ("name", "packaging"),
    [
        ("x.txt.gz", "plain"),
        ("x.zip", "plain"),
        ("x.txt", "gzip"),
        ("a.zip.d/x.txt", "plain"),
        ("b.gz/x.txt", "plain"),
    ],
)
def test_misleading_names(tmp_path, data, expected, name, packaging):
    path = tmp_path / name
    path.parent.mkdir(exist_ok=True)
    path.write_bytes(packed(data, packaging))
    assert_same_data(read_file(path), expected)
    assert_same_data(read_file(str(path)), expected)


def test_path_gzip_packaging_names_the_file(tmp_path):
    content = package(FIXTURE, "gzip-fname", tmp_path).read_bytes()
    # FNAME flag of the gzip header
    assert content[3] & 0x08


def test_header_keeps_line_endings(tmp_path):
    path = tmp_path / "23andme_win.txt.gz"
    path.write_bytes(gzip.compress((INPUT / "23andme_win.txt").read_bytes()))
    with InputSource(path) as source:
        header = source.header()
    assert header.first_line == "# 23andMe\r\n"
    assert header.comments.endswith("# rsid\tchromosome\tposition\tgenotype\r\n")


def test_header_without_data_section_stops_at_end():
    with InputSource(b"[Header]\nContent\t\tinvented\n") as source:
        header = source.header()
    assert header == Header("[Header]\n", "[Header]\nContent\t\tinvented\n")


def test_header_without_data_section_is_parse_error():
    with pytest.raises(ParseError, match="GSA"):
        read_file(b"[Header]\nContent\t\tinvented\n")


class RecordingBytesIO(io.BytesIO):
    """BytesIO that counts the bytes read from it and the furthest offset reached."""

    def __init__(self, data):
        super().__init__(data)
        self.count = 0
        self.furthest = 0

    def read(self, size=-1):
        data = super().read(size)
        self.count += len(data)
        self.furthest = max(self.furthest, self.tell())
        return data

    def seek(self, offset, whence=io.SEEK_SET):
        position = super().seek(offset, whence)
        self.furthest = max(self.furthest, position)
        return position


def invented_rows(count):
    """Return 23andMe rows of invented SNPs from rs101 on."""
    return "".join(f"rs{100 + i}\t1\t{i}\tAA\n" for i in range(1, count + 1)).encode()


@pytest.mark.parametrize("packaging", PACKAGINGS)
def test_header_leaves_body_unread(data, packaging):
    content = packed(data + invented_rows(200_000), packaging)
    stream = RecordingBytesIO(content)
    with InputSource(stream) as source:
        header = source.header()
    assert header.first_line == data.splitlines(keepends=True)[0].decode()
    assert stream.count < len(content) // 2
    # a zip archive is found from its end, but nothing else is looked for there
    if packaging != "zip":
        assert stream.furthest < len(content) // 2


def test_caller_gzip_file_read_once(data):
    backing = RecordingBytesIO(gzip.compress(data + invented_rows(200_000)))
    with gzip.GzipFile(fileobj=backing) as stream:
        s = SNPs(stream)
        assert stream.tell() == 0
    assert s.source == "23andMe"
    assert s.count == 8 + 200_000
    assert backing.count < 1.5 * len(backing.getvalue())


@pytest.mark.parametrize("packaging", PACKAGINGS)
def test_open_gives_independent_streams(data, packaging):
    with InputSource(packed(data, packaging)) as source:
        with source.open() as first, source.open() as second:
            assert first.read(10) == data[:10]
            assert second.read(20) == data[:20]
            assert first.read() == data[10:]
            assert second.read() == data[20:]


@pytest.mark.parametrize("packaging", PACKAGINGS)
def test_file_object_read_from_its_position(data, expected, packaging):
    junk = b"invented junk\x00\x1f\x8bPK\x03\x04"
    content = junk + packed(data, packaging)
    stream = io.BytesIO(content)
    stream.seek(len(junk))
    assert_same_data(read_file(stream), expected)
    assert not stream.closed
    assert stream.tell() == len(junk)
    assert stream.getvalue() == content


AFTER_ERROR = {
    "unrecognized": (lambda data: UNRECOGNIZED, UnsupportedFormatError),
    "several-data-files": (
        lambda data: zip_bytes({"a.txt": data, "b.txt": data}),
        UnsupportedFormatError,
    ),
    "truncated-gzip": (lambda data: gzip.compress(data)[:-20], ParseError),
}


@pytest.mark.parametrize(("make", "error"), AFTER_ERROR.values(), ids=AFTER_ERROR)
def test_file_object_position_restored_after_error(data, make, error):
    junk = b"invented junk"
    stream = io.BytesIO(junk + make(data))
    stream.seek(len(junk))
    with pytest.raises(error):
        read_file(stream)
    assert not stream.closed
    assert stream.tell() == len(junk)


@pytest.mark.parametrize("spill", [False, True], ids=["in-memory", "on-disk"])
@pytest.mark.parametrize("packaging", PACKAGINGS)
def test_non_seekable_stream(monkeypatch, data, expected, packaging, spill):
    if spill:
        monkeypatch.setattr(snps.io.readers._source, "_SPOOL_SIZE", 16)
    stream = NonSeekable(packed(data, packaging))
    assert_same_data(read_file(stream), expected)
    assert not stream.closed


# entries of each archive, None standing for the data file
@pytest.mark.parametrize(
    "entries",
    [
        {"folder/": b"", "folder/a.txt": None},
        {"__MACOSX/._a.txt": b"\x00\x05\x16\x07", "a.txt": None},
        {".DS_Store": b"\x00\x00\x00\x01", "a.txt": None},
        {"folder/.x": b"\x00\x00\x00\x01", "folder/a.txt": None},
        {"readme.txt": UNRECOGNIZED, f"genotypes/{GFG}": None},
    ],
    ids=["directory", "macosx", "ds-store", "nested-dot-file", "genes-for-good"],
)
def test_zip_reads_the_data_file(data, expected, entries):
    archive = zip_bytes({name: data if v is None else v for name, v in entries.items()})
    assert_same_data(read_file(archive), expected)


def test_zip_with_several_data_files(data):
    with pytest.raises(UnsupportedFormatError, match="several") as info:
        read_file(zip_bytes({"a.txt": data, "b.txt": data}))
    assert "'a.txt'" in str(info.value)
    assert "'b.txt'" in str(info.value)


def test_zip_with_two_genes_for_good_files(data):
    entries = {f"one/{GFG}": data, f"two/{GFG}": data}
    with pytest.raises(UnsupportedFormatError, match="several"):
        read_file(zip_bytes(entries))


@pytest.mark.parametrize(
    "entries",
    [{}, {"folder/": b"", "__MACOSX/._a.txt": b"\x00", ".DS_Store": b"\x00"}],
    ids=["empty", "metadata-only"],
)
def test_zip_with_no_data_file(entries):
    with pytest.raises(UnsupportedFormatError, match="no data file"):
        read_file(zip_bytes(entries))


@pytest.mark.parametrize(
    ("fields", "match"),
    [({"flags": 0x1}, "encrypted"), ({"method": 99}, "not supported")],
    ids=["encrypted", "unknown-method"],
)
def test_zip_member_that_cannot_be_read(data, fields, match):
    with pytest.raises(UnsupportedFormatError, match=match):
        read_file(stored_zip(data, **fields))


@pytest.mark.parametrize("packaging", ["zip-stored", "zip-deflated"])
def test_zipped_vcf(tmp_path, packaging):
    vcf = INPUT / "testvcf.vcf"
    path = package(vcf, packaging, tmp_path)
    expected = read_file(vcf)
    assert_same_data(read_file(path), expected)
    assert_same_data(read_file(path.read_bytes()), expected)


@pytest.mark.parametrize(
    ("compression", "make", "cause"), CORRUPT.values(), ids=CORRUPT
)
@pytest.mark.parametrize("kind", ["path", "bytes", "file"])
def test_corrupt_compressed_data(tmp_path, data, kind, compression, make, cause):
    content = make(data)
    path = tmp_path / "corrupt.txt"
    path.write_bytes(content)
    if kind == "path":
        file, name = path, f"'{path}'"
    elif kind == "bytes":
        file, name = content, "bytes input"
    else:
        file, name = io.BytesIO(content), "file object input"
    with pytest.raises(ParseError) as info:
        read_file(file)
    assert_corrupt_data_error(info.value, name, compression, cause)


CALLERS = {
    "read_file": read_file,
    "SNPs": SNPs,
    "SNPs-detect": partial(SNPs, only_detect_source=True),
    "Reader": lambda file: Reader(file).read(),
    "Reader-detect": lambda file: Reader(file, only_detect_source=True).read(),
}


@pytest.mark.parametrize(
    "fault", ["truncated-gzip", "zip-bad-crc", "zip-undecodable-local-name"]
)
@pytest.mark.parametrize("caller", CALLERS.values(), ids=CALLERS)
def test_corrupt_data_error_for_each_caller(data, caller, fault):
    compression, make, cause = CORRUPT[fault]
    with pytest.raises(ParseError) as info:
        caller(make(data))
    assert_corrupt_data_error(info.value, "bytes input", compression, cause)


@pytest.mark.parametrize(
    ("packaging", "cause"), [("gzip", EOFError), ("zip", zipfile.BadZipFile)]
)
@pytest.mark.parametrize(
    "caller", [read_file, SNPs, CALLERS["Reader"]], ids=["read_file", "SNPs", "Reader"]
)
def test_corrupt_data_found_after_header(data, caller, packaging, cause):
    rows = "".join(f"rs{i}\t1\t{i}\tAA\n" for i in range(1, 20_001))
    plain = data + rows.encode()
    if packaging == "gzip":
        content = gzip.compress(plain)[:-20]
    else:
        content = stored_zip(plain, crc_change=0xFFFF)
    with InputSource(content) as source:
        assert source.header().first_line == data.splitlines(keepends=True)[0].decode()
    with pytest.raises(ParseError) as info:
        caller(content)
    assert_corrupt_data_error(info.value, "bytes input", packaging, cause)


def test_gzip_with_trailing_bytes(data, expected):
    assert_same_data(read_file(gzip.compress(data) + b"extra data"), expected)


@pytest.mark.parametrize("padding", [0, 8], ids=["adjacent", "zero-padded"])
def test_gzip_members_are_joined(data, expected, padding):
    middle = data.index(b"rs5")
    content = gzip.compress(data[:middle]) + b"\x00" * padding
    assert_same_data(read_file(content + gzip.compress(data[middle:])), expected)


@pytest.mark.parametrize("padding", range(1, 7))
def test_gzip_padding_split_across_reads(monkeypatch, data, expected, padding):
    monkeypatch.setattr(snps.io.readers._source, "_CHUNK_SIZE", 3)
    middle = data.index(b"rs5")
    content = gzip.compress(data[:middle]) + b"\x00" * padding
    assert_same_data(read_file(content + gzip.compress(data[middle:])), expected)


@pytest.mark.parametrize(
    "trailer", [b"\x00" * 8, b"\x00" * 8 + b"extra data"], ids=["zeros", "zeros-junk"]
)
def test_gzip_trailing_padding(data, expected, trailer):
    assert_same_data(read_file(gzip.compress(data) + trailer), expected)


def count_opens(monkeypatch):
    """Record each stream opened from an input source."""
    opens = []
    open_ = InputSource.open

    def recording_open(self):
        opens.append(self.name)
        return open_(self)

    monkeypatch.setattr(InputSource, "open", recording_open)
    return opens


RETRIED = {
    # the python engine of the last retry rejects the field after the quotes
    "generic": (
        'rsid,chromosome,position,genotype\nrs1,1,101,AA\nrs2,1,102,"A"A\n',
        "rs{i},1,{i},AA\n",
    ),
    "FTDNA": (
        'RSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","101","AA"\n"rs2","1","102","A"A\n',
        '"rs{i}","1","{i}","AA"\n',
    ),
}


@pytest.mark.parametrize(("start", "row"), RETRIED.values(), ids=RETRIED)
def test_corrupt_data_is_not_retried(monkeypatch, start, row):
    text = start + "".join(row.format(i=i) for i in range(1000, 11_000))
    content = gzip.compress(text.encode())[:-20]
    opens = count_opens(monkeypatch)
    with pytest.raises(ParseError) as info:
        read_file(content)
    assert_corrupt_data_error(info.value, "bytes input", "gzip", EOFError)
    # the header, then one parse
    assert len(opens) <= 2


@pytest.mark.parametrize("packaging", PACKAGINGS)
def test_detect_then_read_same_file_object(data, packaging):
    stream = io.BytesIO(packed(data, packaging))
    assert SNPs(stream, only_detect_source=True).source == "23andMe"
    assert_matches_case(SNPs(stream), CASE)


@pytest.mark.parametrize("packaging", PACKAGINGS)
def test_read_twice(data, packaging):
    content = packed(data, packaging)
    assert_same_data(read_file(content), read_file(content))
    stream = io.BytesIO(content)
    assert_same_data(read_file(stream), read_file(stream))


def test_snps_reads_path_and_file_object():
    assert_matches_case(SNPs(FIXTURE), CASE)
    with FIXTURE.open("rb") as f:
        assert_matches_case(SNPs(f), CASE)
        assert not f.closed


def test_snps_gives_empty_object_for_missing_path(tmp_path):
    s = SNPs(tmp_path / "missing.txt")
    assert s.count == 0
    assert s.source == ""
