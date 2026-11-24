"""The input of one read: a raw data file given as a path, bytes or a file object."""

from __future__ import annotations

import errno
import io
import os
import tempfile
import zipfile
import zlib
from collections.abc import Callable, Iterator
from contextlib import ExitStack, contextmanager
from typing import IO, TYPE_CHECKING, BinaryIO, Literal, Self, cast

from snps.exceptions import ParseError, UnsupportedFormatError
from snps.io.readers._common import extract_comments
from snps.io.readers._registry import Header

if TYPE_CHECKING:
    from types import TracebackType

    from _typeshed import WriteableBuffer

_GZIP_MAGIC = b"\x1f\x8b"
# signatures of a zip local file header and of an end of central directory record
_ZIP_MAGIC = (b"PK\x03\x04", b"PK\x05\x06")
_GENES_FOR_GOOD = "GFG_filtered_unphased_genotypes_23andMe.txt"
# general purpose flag bit of an encrypted zip entry
_ENCRYPTED = 0x1
_CHUNK_SIZE = 64 * 1024
_SPOOL_SIZE = 16 * 1024 * 1024
# errors of decompressing corrupt or truncated data
_DECOMPRESSION_ERRORS: tuple[type[Exception], ...] = (
    zipfile.BadZipFile,
    zlib.error,
    EOFError,
)
try:
    import lzma
except ImportError:  # Python built without lzma
    pass
else:
    _DECOMPRESSION_ERRORS += (lzma.LZMAError,)


def can_read(file: object) -> bool:
    """Return whether ``file`` is bytes, a binary file object or an existing path.

    A file object that reads text is not a binary file object.
    """
    if isinstance(file, (str, os.PathLike)):
        return os.path.exists(file)
    if isinstance(file, bytes):
        return True
    if not callable(getattr(file, "read", None)):
        return False
    return bool(getattr(file, "closed", False)) or not _reads_text(file)


def _reads_text(file: object) -> bool:
    """Return whether an open file object reads text, by reading nothing from it."""
    if isinstance(file, io.TextIOBase):
        return True
    read: Callable[[int], object] = getattr(file, "read")
    return isinstance(read(0), str)


class InputSource:
    """The data of one raw data file, which may be gzip or zip compressed.

    Parameters
    ----------
    file : str, os.PathLike, bytes or binary file object
        file to read, as for ``SNPs``

    Attributes
    ----------
    name : str
        description of the input for messages
    compression : {"gzip", "zip"} or None
        compression of the data, found from its leading bytes

    Raises
    ------
    TypeError
        if ``file`` is none of those kinds, or is a text stream
    ValueError
        if ``file`` is a closed file object
    FileNotFoundError
        if no file exists at the path
    IsADirectoryError
        if the path is a directory
    UnsupportedFormatError
        if a zip archive does not hold exactly one data file to read, or holds it
        encrypted or compressed in a way that is not supported
    ParseError
        if the compressed data is corrupt or truncated
    """

    compression: Literal["gzip", "zip"] | None

    def __init__(self, file: str | os.PathLike[str] | bytes | BinaryIO) -> None:
        self._member: zipfile.ZipInfo | None = None
        self._errors: list[ParseError] = []
        with ExitStack() as resources:
            self.name, self._open_raw = _opener(file, resources)
            self.compression = self._detect()
            self._resources = resources.pop_all()

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        self.close()

    def close(self) -> None:
        """Release what the source holds, and return a file object to its position."""
        self._resources.close()

    def raised(self, exc: BaseException) -> bool:
        """Return whether ``exc`` is an error the source raised for corrupt data."""
        return any(exc is error for error in self._errors)

    def header(self) -> Header:
        """Read the first line and leading comments of the data.

        Returns
        -------
        Header
            start of the data, decoded as UTF-8 with its line endings kept

        Raises
        ------
        UnicodeDecodeError
            if the start of the data is not UTF-8 text
        ParseError
            if the compressed data is corrupt or truncated
        """
        with self.open() as stream:
            first_line, comments, _ = extract_comments(stream)
        return Header(first_line, comments)

    def open(self) -> IO[bytes]:
        """Open a new stream of the plain data, at its start, for the caller to close.

        Returns
        -------
        binary file object
            plain (decompressed) data

        Raises
        ------
        ParseError
            if the compressed data is corrupt or truncated; also raised while reading
        """
        if self.compression is None:
            return self._open_raw()
        with ExitStack() as resources:
            raw = resources.enter_context(self._open_raw())
            read: Callable[[int], bytes]
            if self._member is None:
                read = _GzipMembers(raw).read
            else:
                with _on_corrupt_zip(self._corrupt):
                    archive = resources.enter_context(zipfile.ZipFile(raw))
                    read = resources.enter_context(archive.open(self._member)).read
            plain = _PlainData(read, resources.pop_all(), self._corrupt)
        return io.BufferedReader(plain, _CHUNK_SIZE)

    def _detect(self) -> Literal["gzip", "zip"] | None:
        """Find the compression of the data, and the zip entry to read."""
        with self._open_raw() as raw:
            lead = raw.read(len(_ZIP_MAGIC[0]))
            if lead.startswith(_GZIP_MAGIC):
                return "gzip"
            if lead not in _ZIP_MAGIC:
                return None
            raw.seek(0)
            try:
                with _on_corrupt_zip(lambda: self._corrupt("zip")):
                    archive = zipfile.ZipFile(raw)
            except NotImplementedError as exc:
                raise UnsupportedFormatError(
                    f"{self.name} is a zip archive of a version that is not supported."
                ) from exc
            with archive:
                self._member = self._find_member(archive)
            return "zip"

    def _find_member(self, archive: zipfile.ZipFile) -> zipfile.ZipInfo:
        """Return the one data file of a zip archive, checking that it can be read."""
        entries = [info for info in archive.infolist() if _is_data_file(info)]
        if not entries:
            raise UnsupportedFormatError(
                f"{self.name} is a zip archive with no data file."
            )
        if len(entries) > 1:
            matches = [info for info in entries if _GENES_FOR_GOOD in info.filename]
            if len(matches) != 1:
                names = ", ".join(repr(info.filename) for info in entries)
                raise UnsupportedFormatError(
                    f"{self.name} is a zip archive with several data files ({names}); "
                    "it must hold one raw data file."
                )
            entries = matches
        member = entries[0]
        if member.flag_bits & _ENCRYPTED:
            raise UnsupportedFormatError(
                f"{self.name} holds {member.filename!r} encrypted; "
                "decrypt it before reading."
            )
        if member.header_offset < 0:
            raise self._corrupt("zip")
        try:
            with _on_corrupt_zip(lambda: self._corrupt("zip")), archive.open(member):
                pass
        except (NotImplementedError, RuntimeError) as exc:
            raise UnsupportedFormatError(
                f"{self.name} holds {member.filename!r} compressed in a way that is "
                f"not supported ({exc})."
            ) from exc
        return member

    def _corrupt(self, compression: str | None = None) -> ParseError:
        """Return the error for compressed data that is corrupt or truncated."""
        error = ParseError(
            f"Cannot read {self.name}: the {compression or self.compression} "
            "compressed data is corrupt or truncated."
        )
        self._errors.append(error)
        return error


def _opener(file: object, resources: ExitStack) -> tuple[str, Callable[[], IO[bytes]]]:
    """Return the name of an input, and a function opening its raw data from the start.

    A file object is returned to its position, or copied if it cannot seek, by
    closing ``resources``.
    """
    if isinstance(file, bytes):
        return "bytes input", lambda: io.BytesIO(file)
    if isinstance(file, (str, os.PathLike)):
        path = os.fspath(file)
        if not os.path.exists(path):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path)
        if os.path.isdir(path):
            raise IsADirectoryError(
                errno.EISDIR, "Is a directory; pass the path of one raw data file", path
            )
        return f"'{path}'", lambda: open(path, "rb")
    if not callable(getattr(file, "read", None)):
        raise TypeError(
            "file must be a path, bytes or a binary file object, "
            f"not {type(file).__name__}"
        )
    stream = cast("IO[bytes]", file)
    if getattr(stream, "closed", False):
        raise ValueError("I/O operation on closed file.")
    if _reads_text(stream):
        raise TypeError("file must be opened in binary mode, not text mode")
    name = getattr(stream, "name", None)
    label = f"'{name}'" if isinstance(name, str) else "file object input"
    seekable = getattr(stream, "seekable", None)
    if callable(seekable) and seekable():
        start = stream.tell()
        resources.callback(stream.seek, start)
    else:
        spool = resources.enter_context(
            tempfile.SpooledTemporaryFile(max_size=_SPOOL_SIZE)
        )
        while chunk := _binary(stream.read(_CHUNK_SIZE)):
            spool.write(chunk)
        stream, start = spool, 0
    return label, lambda: io.BufferedReader(_Window(stream, start), _CHUNK_SIZE)


def _binary(data: object) -> bytes:
    """Return data read from a file object, which must be bytes."""
    if not isinstance(data, bytes):
        raise TypeError(
            "file must be a binary file object, but its read returned "
            f"{type(data).__name__}"
        )
    return data


@contextmanager
def _on_corrupt_zip(error: Callable[[], ParseError]) -> Iterator[None]:
    """Raise ``error()`` from an error opening a zip archive or entry that is corrupt.

    BadZipFile and ValueError mean corrupt data, and so does an OSError for an
    invalid argument: a corrupt offset fails to seek, with ValueError in bytes and
    EINVAL in a file. Any other OSError propagates.
    """
    try:
        yield
    except (zipfile.BadZipFile, ValueError, OSError) as exc:
        if isinstance(exc, OSError) and exc.errno != errno.EINVAL:
            raise
        raise error() from exc


def _is_data_file(info: zipfile.ZipInfo) -> bool:
    """Return whether a zip entry may be the data file, not a directory or metadata."""
    name = info.filename
    return not (
        info.is_dir()
        or name.startswith("__MACOSX/")
        or name.rsplit("/", 1)[-1].startswith(".")
    )


class _Window(io.RawIOBase):
    """The data of a seekable binary file from an offset to its end.

    Each window keeps its own position, and closing it leaves the file open.
    """

    def __init__(self, file: IO[bytes], start: int) -> None:
        super().__init__()
        self._file = file
        self._start = start
        self._position = 0

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self._position

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        if whence == io.SEEK_SET:
            position = offset
        elif whence == io.SEEK_CUR:
            position = self._position + offset
        elif whence == io.SEEK_END:
            position = self._file.seek(0, io.SEEK_END) - self._start + offset
        else:
            raise ValueError(f"invalid whence ({whence})")
        if position < 0:
            raise OSError(errno.EINVAL, "negative seek position", position)
        self._position = position
        return position

    def readinto(self, buffer: WriteableBuffer) -> int:
        self._file.seek(self._start + self._position)
        with memoryview(buffer) as view, view.cast("B") as target:
            data = _binary(self._file.read(len(target)))
            target[: len(data)] = data
        self._position += len(data)
        return len(data)


class _PlainData(io.RawIOBase):
    """Plain data read from compressed data; closing it closes what it reads from."""

    def __init__(
        self,
        read: Callable[[int], bytes],
        resources: ExitStack,
        corrupt: Callable[[], ParseError],
    ) -> None:
        super().__init__()
        self._read = read
        self._resources = resources
        self._corrupt = corrupt

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: WriteableBuffer) -> int:
        with memoryview(buffer) as view, view.cast("B") as target:
            try:
                data = self._read(len(target))
            except _DECOMPRESSION_ERRORS as exc:
                raise self._corrupt() from exc
            except OSError as exc:
                # bz2 reports corrupt data as an OSError without an errno
                if exc.errno is not None:
                    raise
                raise self._corrupt() from exc
            target[: len(data)] = data
        return len(data)

    def close(self) -> None:
        if not self.closed:
            try:
                self._resources.close()
            finally:
                super().close()


class _GzipMembers:
    """Reader of the plain data of the gzip members at the start of a stream.

    Zero bytes after a member are padding, as the gzip module reads them; other
    bytes that do not start another member end the data.
    """

    def __init__(self, stream: IO[bytes]) -> None:
        self._stream = stream
        self._input = b""
        self._decompressor = zlib.decompressobj(wbits=31)
        self._in_member = False
        self._done = False

    def read(self, size: int) -> bytes:
        """Return up to ``size`` bytes of plain data, or no bytes at the end.

        Raises
        ------
        zlib.error
            if a member is corrupt
        EOFError
            if the data ends inside a member
        """
        while size > 0 and not self._done:
            if not self._in_member and not self._next_member():
                self._done = True
                break
            if not self._input:
                self._input = self._stream.read(_CHUNK_SIZE)
                if not self._input:
                    raise EOFError("gzip data ended inside a member")
            data = self._decompressor.decompress(self._input, size)
            if self._decompressor.eof:
                self._input = self._decompressor.unused_data
                self._in_member = False
            else:
                self._input = self._decompressor.unconsumed_tail
            if data:
                return data
        return b""

    def _next_member(self) -> bool:
        """Start decompressing the next member, returning whether one follows."""
        self._input = self._input.lstrip(b"\x00")
        while len(self._input) < len(_GZIP_MAGIC):
            chunk = self._stream.read(_CHUNK_SIZE)
            if not chunk:
                break
            self._input = (self._input + chunk).lstrip(b"\x00")
        if not self._input.startswith(_GZIP_MAGIC):
            return False
        self._decompressor = zlib.decompressobj(wbits=31)
        self._in_member = True
        return True
