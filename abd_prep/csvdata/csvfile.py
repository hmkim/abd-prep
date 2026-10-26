"""CSV reading: delimiter detection, header checks, row streaming (F19).

Rows are streamed with :mod:`csv` so a file close to the 100 MB ABD limit
never has to be loaded whole. Line numbers come from ``reader.line_num``,
which counts physical lines and therefore stays right across quoted
fields that span several lines.

Finding codes:

``CSV_NO_HEADER``            error    empty file / no header row
``CSV_HEADER_EMPTY_NAME``    error    a header cell is blank
``CSV_HEADER_DUPLICATE``     error    two header cells have the same name
``CSV_ROW_LENGTH``           error    a row has more/fewer cells than the header (row excluded from clean output)
``CSV_DELIMITER``            info     delimiter is not ',' (clean writes ',')
``CSV_DELIMITER_GUESSED``    warning  Sniffer could not decide; ',' assumed -- pass --delimiter
``CSV_BLANK_ROWS``           info     completely empty rows skipped
"""
import csv
import io
import os
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Tuple

from ..common.result import InputError, ValidationResult

SNIFF_DELIMITERS = ",;\t|"
SNIFF_BYTES = 64 * 1024
LIST_LIMIT = 10

# The csv module refuses fields above 128 KB by default; a long free-text
# cell is legitimate, so allow anything up to the ABD file limit.
csv.field_size_limit(100 * 1024 * 1024)

_DELIMITER_NAMES: Dict[str, str] = {
    "tab": "\t", "\\t": "\t", "t": "\t",
    "comma": ",", "semicolon": ";", "pipe": "|",
}


def _csv_failure(path: str, lineno: int, exc: csv.Error) -> InputError:
    return InputError(f"{os.path.basename(path)}: malformed CSV near line {lineno}: {exc}. "
                      "Check quoting (a stray '\"' usually) or pass --delimiter")


def parse_delimiter(text: Optional[str]) -> Optional[str]:
    """``--delimiter`` value -> one character (accepts 'tab', '\\t', ';', ',', '|')."""
    if text is None:
        return None
    key = text.strip().lower()
    if key in _DELIMITER_NAMES:
        return _DELIMITER_NAMES[key]
    if len(text) == 1:
        return text
    raise InputError(f"--delimiter must be a single character or one of tab/comma/semicolon/pipe, got {text!r}")


def sniff_delimiter(sample: str) -> Optional[str]:
    """Return the delimiter detected in ``sample``, or None if undecidable.

    :class:`csv.Sniffer` first. It gives up on short files and on files
    with an irregular row, so when it fails the header line is examined:
    if exactly one candidate delimiter occurs there, that is the answer
    (header cells rarely contain stray delimiters). Otherwise None and
    the caller assumes ',' with a warning.
    """
    if not sample.strip():
        return None
    try:
        return csv.Sniffer().sniff(sample, delimiters=SNIFF_DELIMITERS).delimiter
    except csv.Error:
        pass
    header = next((ln for ln in sample.splitlines() if ln.strip()), "")
    present = [d for d in SNIFF_DELIMITERS if d in header]
    if len(present) == 1:
        return present[0]
    return None


@dataclass
class CsvLayout:
    delimiter: str
    header: List[str]
    header_lineno: int = 1
    sniffed: bool = False


@dataclass
class RowStats:
    rows: int = 0               # data rows with the right length
    bad_length: List[Tuple[int, int]] = field(default_factory=list)   # (lineno, n_cells)
    blank: int = 0


def read_layout(path: str, codec: str, delimiter: Optional[str], result: ValidationResult) -> Optional[CsvLayout]:
    """Decide the delimiter and read/check the header row. None when there is no header."""
    with open(path, "r", encoding=codec, newline="") as fh:
        sample = fh.read(SNIFF_BYTES)
    if not sample.strip():
        result.error("CSV_NO_HEADER", "file is empty; a header row is required", "file")
        return None
    sniffed = False
    if delimiter is None:
        guess = sniff_delimiter(sample)
        if guess is None:
            delimiter = ","
            result.warning("CSV_DELIMITER_GUESSED",
                           "could not detect the delimiter; ',' assumed. Pass --delimiter ';' (or tab) if the "
                           "file uses another one", "file")
        else:
            delimiter, sniffed = guess, True
    if delimiter != ",":
        name = "tab" if delimiter == "\t" else repr(delimiter)
        how = "detected" if sniffed else "given"
        result.info("CSV_DELIMITER", f"delimiter {name} ({how}); clean writes ','", "file")

    reader = csv.reader(io.StringIO(sample), delimiter=delimiter)
    header: List[str] = []
    try:
        for row in reader:
            if any(c.strip() for c in row):
                header = row
                break
    except csv.Error as exc:
        raise _csv_failure(path, reader.line_num, exc)
    if not header:
        result.error("CSV_NO_HEADER", "no header row found", "file")
        return None
    layout = CsvLayout(delimiter=delimiter, header=header, header_lineno=reader.line_num, sniffed=sniffed)
    names = [h.strip() for h in header]
    empties = [i + 1 for i, h in enumerate(names) if not h]
    if empties:
        result.error("CSV_HEADER_EMPTY_NAME",
                     f"header column(s) {', '.join(str(i) for i in empties)} have no name", "header")
    seen: Dict[str, int] = {}
    for h in names:
        seen[h] = seen.get(h, 0) + 1
    dups = [h for h, n in seen.items() if n > 1 and h]
    if dups:
        result.error("CSV_HEADER_DUPLICATE", f"duplicate header name(s): {', '.join(repr(d) for d in dups)}", "header")
    return layout


def iter_rows(path: str, codec: str, layout: CsvLayout, stats: RowStats) -> Iterator[Tuple[int, List[str]]]:
    """Yield ``(lineno, cells)`` for every data row with the header's width.

    Rows of another width are counted in ``stats.bad_length`` and skipped;
    blank rows are counted and skipped.
    """
    width = len(layout.header)
    with open(path, "r", encoding=codec, newline="") as fh:
        reader = csv.reader(fh, delimiter=layout.delimiter)
        header_seen = False
        try:
            for row in reader:
                if not header_seen:
                    if any(c.strip() for c in row):
                        header_seen = True
                    continue
                if not row or not any(c.strip() for c in row):
                    stats.blank += 1
                    continue
                if len(row) != width:
                    stats.bad_length.append((reader.line_num, len(row)))
                    continue
                stats.rows += 1
                yield reader.line_num, row
        except csv.Error as exc:
            raise _csv_failure(path, reader.line_num, exc)


def report_row_stats(stats: RowStats, width: int, result: ValidationResult) -> None:
    if stats.bad_length:
        shown = ", ".join(f"line {ln} ({n} cells)" for ln, n in stats.bad_length[:LIST_LIMIT])
        more = f" ... ({len(stats.bad_length)} rows)" if len(stats.bad_length) > LIST_LIMIT else ""
        result.error("CSV_ROW_LENGTH",
                     f"{len(stats.bad_length)} row(s) do not have the header's {width} cells: {shown}{more}; "
                     "an unquoted delimiter inside a value is the usual cause. These rows are not written by clean",
                     "rows")
    if stats.blank:
        result.info("CSV_BLANK_ROWS", f"{stats.blank} blank row(s) skipped", "rows")
