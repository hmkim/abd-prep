"""File hashing and the two ABD upload constraints shared by every group.

ABD rejects uploads of 100 MB or more and files that are not UTF-8
(plan §8.1). Both checks stream the file so a large CSV never has to be
loaded into memory.
"""
import codecs
import hashlib
import os
from typing import Optional, Tuple

#: ABD upload limit: files must be smaller than this (plan §8.1).
ABD_MAX_BYTES = 100 * 1024 * 1024

_CHUNK = 1024 * 1024


def sha256_file(path: str) -> str:
    """Hex SHA-256 of a file, streamed in 1 MiB chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def file_size(path: str) -> int:
    return os.path.getsize(path)


def exceeds_abd_limit(path: str) -> bool:
    return file_size(path) >= ABD_MAX_BYTES


def check_utf8(path: str, encoding: str = "utf-8") -> Tuple[bool, Optional[str]]:
    """Return ``(ok, detail)`` after streaming the file through a decoder.

    ``detail`` carries the decoder's error message (with byte offset) when
    the file is not valid in the given encoding. A UTF-8 BOM is accepted
    (``utf-8-sig`` semantics) because Excel writes one; callers strip it.
    """
    codec = "utf-8-sig" if encoding.lower().replace("_", "-") in ("utf-8", "utf8") else encoding
    try:
        decoder = codecs.getincrementaldecoder(codec)(errors="strict")
    except LookupError:
        return False, f"unknown encoding {encoding!r}"
    offset = 0            # bytes handed to the decoder *before* the current chunk
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(_CHUNK), b""):
                decoder.decode(chunk)
                offset += len(chunk)
            decoder.decode(b"", final=True)
    except UnicodeDecodeError as exc:
        # exc.start is relative to the decoder's pending bytes + this chunk;
        # offset has not been advanced for the failing chunk yet.
        return False, f"not valid {encoding}: {exc.reason} near byte {offset + exc.start}"
    return True, None


def has_utf8_bom(path: str) -> bool:
    with open(path, "rb") as fh:
        return fh.read(3) == codecs.BOM_UTF8
