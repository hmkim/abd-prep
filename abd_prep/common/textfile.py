"""Pre-checks every text input (FASTA, CSV) goes through before parsing.

The two ABD upload constraints (plan §8.1) -- a file must be smaller
than 100 MB and encoded as UTF-8 -- plus the extension and encoding
conventions of plan §3.1 live here so the ``fasta`` and ``csv`` groups
report them identically.

Finding codes (all ``where="file"``):

``FILE_TOO_LARGE``          error    >= 100 MB (ABD rejects the upload)
``FILE_NOT_UTF8``           error    not valid UTF-8 and no ``--encoding`` given
``FILE_ENCODING``           error    not valid in the ``--encoding`` the user named
``FILE_ENCODING_INPUT``     warning  (validate) input is in another encoding; ABD needs UTF-8
``FILE_ENCODING_CONVERTED`` info     (clean) input decoded as ``--encoding``, output written as UTF-8
``FILE_ENCODING_ASCII``     info     ``--encoding`` given but the file is plain ASCII (already UTF-8)
``FILE_BOM``                info     UTF-8 byte-order mark present (removed by ``clean``)
``FILE_EXTENSION``          warning  extension unusual for the group; the file is still parsed

An Excel workbook (``.xlsx``/``.xls``/``.xlsm``) is an input error (exit 2):
the standard library cannot read it, and the message tells the user to
save it as "CSV UTF-8" from Excel.
"""
import codecs
import os
from typing import Optional, Sequence

from .hashing import ABD_MAX_BYTES, check_utf8, exceeds_abd_limit, has_utf8_bom
from .result import InputError, ValidationResult

UTF8_NAMES = ("utf-8", "utf8", "utf_8", "utf-8-sig", "utf8-sig", "utf_8_sig")
EXCEL_EXTENSIONS = (".xlsx", ".xls", ".xlsm", ".xlsb")

FASTA_EXTENSIONS = (".fasta", ".fa", ".faa", ".fas")
CSV_EXTENSIONS = (".csv", ".tsv", ".txt")


def is_utf8_name(encoding: str) -> bool:
    return encoding.strip().lower() in UTF8_NAMES


def normalize_encoding(encoding: Optional[str]) -> str:
    """Return the codec name to open the file with, or raise ``InputError``.

    UTF-8 (the default) becomes ``utf-8-sig`` so a BOM written by Excel is
    stripped transparently. Any other name must be a codec Python knows.
    """
    enc = (encoding or "utf-8").strip()
    if is_utf8_name(enc):
        return "utf-8-sig"
    try:
        codecs.lookup(enc)
    except LookupError:
        raise InputError(f"unknown encoding {enc!r} (examples: utf-8, cp949, euc-kr, latin-1)")
    return enc


def check_extension(path: str, allowed: Sequence[str], group: str, result: ValidationResult) -> None:
    """Warn about an unusual extension; refuse Excel workbooks outright."""
    lower = path.lower()
    if lower.endswith(EXCEL_EXTENSIONS):
        raise InputError(f"{os.path.basename(path)!r} is an Excel workbook, which this tool cannot read. "
                         "In Excel use File > Save As > 'CSV UTF-8 (Comma delimited)' and run the csv group "
                         "on that file.")
    if not lower.endswith(tuple(allowed)):
        result.warning("FILE_EXTENSION",
                       f"{os.path.basename(path)!r} does not have a usual {group} extension "
                       f"({', '.join(allowed)}); the file is still parsed as {group}", "file")


def precheck(path: str, encoding: Optional[str], result: ValidationResult, action: str = "validate") -> Optional[str]:
    """Run the size/encoding checks and return the codec to decode with.

    Returns ``None`` when the file cannot be decoded (an error finding has
    been recorded and the caller must stop): the input is not UTF-8 and no
    ``--encoding`` was given, or it is not valid in the encoding given.
    ``action`` is ``"clean"`` or ``"validate"`` and only changes the wording
    of the encoding finding.
    """
    if exceeds_abd_limit(path):
        result.error("FILE_TOO_LARGE",
                     f"{os.path.getsize(path)} bytes; ABD rejects uploads of {ABD_MAX_BYTES // (1024 * 1024)} MB "
                     "or more", "file")
    codec = normalize_encoding(encoding)
    if codec == "utf-8-sig":
        ok, detail = check_utf8(path)
        if not ok:
            result.error("FILE_NOT_UTF8",
                         f"{detail}; ABD requires UTF-8. If the file is in another encoding, pass it with "
                         "--encoding (e.g. --encoding cp949) and run clean to convert", "file")
            return None
        if has_utf8_bom(path):
            if action == "clean":
                result.info("FILE_BOM", "UTF-8 byte-order mark present; removed in the output", "file")
            else:
                result.info("FILE_BOM", "UTF-8 byte-order mark present (Excel writes one); clean removes it", "file")
        return codec
    ok, detail = check_utf8(path, codec)
    if not ok:
        result.error("FILE_ENCODING", f"{detail}; check the --encoding value", "file")
        return None
    ascii_only, _ = check_utf8(path, "ascii")
    if ascii_only:
        # identical bytes in every ASCII-compatible encoding: nothing to convert, nothing to warn about
        result.info("FILE_ENCODING_ASCII", f"--encoding {codec} given but the file is plain ASCII, which is already "
                    "valid UTF-8", "file")
        return "utf-8-sig"
    if action == "clean":
        result.info("FILE_ENCODING_CONVERTED", f"input decoded as {codec}; output is written as UTF-8", "file")
    else:
        result.warning("FILE_ENCODING_INPUT",
                       f"file decodes as {codec} but ABD requires UTF-8; run clean with --encoding {codec} "
                       "to convert", "file")
    return codec
