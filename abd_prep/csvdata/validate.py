"""CSV validation rules F20-F21 (plan §3.3.3), the ``clean`` writer and ``inspect``.

One pass over the rows does both validation and, when an output path is
given, writing. Values are copied unchanged; only the *shape* changes:
UTF-8 without BOM, ``,`` delimiter, minimal quoting, LF line ends, the
mapped standard columns first (renamed to their standard names) and the
unmapped source columns after them (dropped with ``--drop-unmapped``).

Finding codes:

``CSV_ID_MISSING``          error    no ``id`` column after mapping
``CSV_SEQUENCE_MISSING``    warning  neither ``sequence`` nor ``heavy``+``light``
``CSV_MISSING_VALUES``      warning  missing tokens, per standard column
``CSV_VALUE_NOT_NUMERIC``   warning  ``value`` cells that are text or censored (``<0.1``)
``CSV_DECIMAL_COMMA``       warning  ``value`` cells like ``1,2``
``CSV_UNIT_MIXED``          warning  more than one distinct ``unit``
``CSV_ID_DUPLICATE``        warning  repeated ``id``
``CSV_SEQUENCE_ALPHABET``   error    sequence cells with letters outside the 20 amino acids
``CSV_SEQUENCE_NUCLEIC``    error    sequence cells that look like DNA/RNA
``CSV_SEQUENCE_LOWERCASE``  info     lower-case sequence cells (left unchanged)
``CSV_ROWS_DROPPED``        info     rows removed by --drop-missing
``CSV_NO_DATA_ROWS``        warning  header only
"""
import csv
import os
from collections import Counter, OrderedDict
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from ..common import residues as R
from ..common.result import ValidationResult
from . import schema as S
from .csvfile import CsvLayout, RowStats, iter_rows, read_layout, report_row_stats

LIST_LIMIT = 8

#: Called for every accepted row with ``{standard_column: cell_text}``.
RowHook = Callable[[int, Dict[str, str]], None]


@dataclass
class CsvOptions:
    mapping: Dict[str, str] = field(default_factory=dict)     # standard -> source header
    delimiter: Optional[str] = None
    drop_missing: bool = False
    drop_unmapped: bool = False


class _Agg:
    """Collect occurrences per key and emit one finding each (keeps reports short)."""

    def __init__(self) -> None:
        self.count: Counter = Counter()
        self.first: Dict[str, List[str]] = {}

    def add(self, key: str, where: str) -> None:
        self.count[key] += 1
        lst = self.first.setdefault(key, [])
        if len(lst) < LIST_LIMIT:
            lst.append(where)

    def describe(self, key: str) -> str:
        n = self.count[key]
        shown = ", ".join(self.first[key])
        return f"{n} cell(s): {shown}" + (" ..." if n > LIST_LIMIT else "")


def process_csv(path: str, codec: str, opts: CsvOptions, result: ValidationResult,
                out_path: Optional[str] = None, row_hook: Optional[RowHook] = None) -> Dict[str, object]:
    """Validate ``path`` and, when ``out_path`` is set, write the cleaned CSV.

    Returns the ``details`` dictionary for report.json.
    """
    layout = read_layout(path, codec, opts.delimiter, result)
    details: Dict[str, object] = {"delimiter": None, "columns": {}, "rows_read": 0, "rows_written": 0,
                                  "output_written": False}
    if layout is None:
        return details
    details["delimiter"] = "tab" if layout.delimiter == "\t" else layout.delimiter
    resolved, unmapped = S.resolve_columns(layout.header, opts.mapping)
    names = [h.strip() for h in layout.header]
    details["columns"] = OrderedDict((std, names[i]) for std, i in sorted(resolved.items(), key=lambda kv: kv[1]))
    details["unmapped"] = [names[i] for i in unmapped]

    for col in S.REQUIRED_COLUMNS:
        if col not in resolved:
            result.error("CSV_ID_MISSING",
                         f"no {col!r} column after mapping; rows cannot be identified. Map one with --map "
                         f"{col}=<header> (header is: {', '.join(repr(n) for n in names)})", "header")
    if not S.has_sequence_columns(resolved):
        result.warning("CSV_SEQUENCE_MISSING",
                       "no 'sequence' (or 'heavy'+'light') column after mapping; cross-referencing with a seed "
                       "FASTA will match by id only", "header")

    std_order = [c for c in S.STANDARD_COLUMNS if c in resolved]
    out_header = std_order + ([] if opts.drop_unmapped else [names[i] for i in unmapped])
    out_indices = [resolved[c] for c in std_order] + ([] if opts.drop_unmapped else unmapped)

    missing = _Agg()
    not_numeric = _Agg()
    decimal_comma = _Agg()
    seq_bad = _Agg()
    seq_nucleic = _Agg()
    seq_lower = _Agg()
    units: Counter = Counter()
    ids: Counter = Counter()
    id_first: Dict[str, int] = {}
    dropped = 0
    stats = RowStats()

    writer = None
    out_fh = None
    rows_written = 0
    if out_path is not None:
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        out_fh = open(out_path, "w", encoding="utf-8", newline="")
        writer = csv.writer(out_fh, delimiter=",", quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
        writer.writerow(out_header)
    try:
        for lineno, row in iter_rows(path, codec, layout, stats):
            row_missing = False
            values: Dict[str, str] = {}
            for std, idx in resolved.items():
                cell = row[idx]
                values[std] = cell
                if S.is_missing(cell):
                    row_missing = True
                    missing.add(std, f"line {lineno}")
                    continue
                if std in S.NUMERIC_COLUMNS:
                    kind, _ = S.classify_number(cell)
                    if kind == S.KIND_DECIMAL_COMMA:
                        decimal_comma.add(std, f"line {lineno} {cell.strip()!r}")
                    elif kind in (S.KIND_CENSORED, S.KIND_TEXT):
                        not_numeric.add(std, f"line {lineno} {cell.strip()!r}")
                elif std == "unit":
                    units[cell.strip()] += 1
                elif std == "id":
                    key = cell.strip()
                    ids[key] += 1
                    id_first.setdefault(key, lineno)
                elif std in S.SEQUENCE_COLUMNS:
                    text = cell.strip()
                    if any(ch.islower() for ch in text):
                        seq_lower.add(std, f"line {lineno}")
                    bad, nucleic = R.check_protein_alphabet(text.upper())
                    if nucleic and len(text) >= 20:
                        seq_nucleic.add(std, f"line {lineno}")
                    elif bad:
                        shown = "".join(sorted({ch for _, ch in bad}))
                        seq_bad.add(std, f"line {lineno} ({shown!r})")
            if opts.drop_missing and row_missing:
                dropped += 1
                continue
            if row_hook is not None:
                row_hook(lineno, values)
            if writer is not None:
                writer.writerow([row[i] for i in out_indices])
                rows_written += 1
    except BaseException:
        if out_fh is not None:
            out_fh.close()
            try:
                os.remove(out_path)          # never leave a half-written output behind
            except OSError:
                pass
        raise
    finally:
        if out_fh is not None and not out_fh.closed:
            out_fh.close()
    details["rows_written"] = rows_written
    details["output_written"] = out_fh is not None

    details["rows_read"] = stats.rows
    report_row_stats(stats, len(layout.header), result)
    if stats.rows == 0:
        result.warning("CSV_NO_DATA_ROWS", "header only; no data rows", "rows")

    for std in sorted(missing.count):
        result.warning("CSV_MISSING_VALUES",
                       f"missing values in {std!r}: {missing.describe(std)}"
                       + ("" if opts.drop_missing else "; --drop-missing removes these rows"), std)
    for std in sorted(decimal_comma.count):
        result.warning("CSV_DECIMAL_COMMA",
                       f"{decimal_comma.describe(std)} contain a ',' inside the number: a European decimal mark "
                       "(1,2 = 1.2) or a thousands separator (1,200)? Save the file with plain '.' decimals and no "
                       "grouping -- values are not converted", std)
    for std in sorted(not_numeric.count):
        result.warning("CSV_VALUE_NOT_NUMERIC",
                       f"{not_numeric.describe(std)} are not plain numbers (assay limits like '<0.1' or text); "
                       "values are kept as written", std)
    if len(units) > 1:
        desc = ", ".join(f"{u!r} x{n}" for u, n in units.most_common())
        result.warning("CSV_UNIT_MIXED", f"{len(units)} different units: {desc}; values are not converted", "unit")
    dups = [(k, n) for k, n in ids.items() if n > 1]
    if dups:
        shown = ", ".join(f"{k!r} x{n}" for k, n in dups[:LIST_LIMIT]) + (" ..." if len(dups) > LIST_LIMIT else "")
        result.warning("CSV_ID_DUPLICATE", f"{len(dups)} id value(s) repeat: {shown}", "id")
    for std in sorted(seq_nucleic.count):
        result.error("CSV_SEQUENCE_NUCLEIC",
                     f"{seq_nucleic.describe(std)} look like DNA/RNA (A/C/G/T/U/N only); amino-acid sequences "
                     "are required", std)
    for std in sorted(seq_bad.count):
        result.error("CSV_SEQUENCE_ALPHABET",
                     f"{seq_bad.describe(std)} contain letters outside the 20 standard amino acids", std)
    for std in sorted(seq_lower.count):
        result.info("CSV_SEQUENCE_LOWERCASE", f"{seq_lower.describe(std)} contain lower-case letters (left as-is)", std)
    if dropped:
        result.info("CSV_ROWS_DROPPED", f"{dropped} row(s) with missing values not written (--drop-missing)", "rows")

    details["missing"] = dict(missing.count)
    details["units"] = dict(units)
    details["ids"] = len(ids)
    details["dropped"] = dropped
    return details


# -- inspect -----------------------------------------------------------------------

def inspect_csv(path: str, codec: str, delimiter: Optional[str], result: ValidationResult,
                mapping: Optional[Dict[str, str]] = None, sample: int = 3) -> Dict[str, object]:
    layout = read_layout(path, codec, delimiter, result)
    if layout is None:
        return {"path": path, "columns": [], "rows": 0}
    names = [h.strip() for h in layout.header]
    resolved, _ = S.resolve_columns(layout.header, mapping or {})
    std_of = {idx: std for std, idx in resolved.items()}
    filled = [0] * len(names)
    samples: List[List[str]] = [[] for _ in names]
    stats = RowStats()
    for _, row in iter_rows(path, codec, layout, stats):
        for i, cell in enumerate(row):
            if not S.is_missing(cell):
                filled[i] += 1
                if len(samples[i]) < sample and cell.strip() not in samples[i]:
                    samples[i].append(cell.strip())
    report_row_stats(stats, len(names), result)
    cols = [{"name": n, "standard": std_of.get(i), "filled": filled[i], "samples": samples[i]}
            for i, n in enumerate(names)]
    return {"path": path, "delimiter": "tab" if layout.delimiter == "\t" else layout.delimiter,
            "rows": stats.rows, "bad_length_rows": len(stats.bad_length), "columns": cols}


def format_inspect(info: Dict[str, object]) -> str:
    lines = [f"== {info['path']}  delimiter: {info.get('delimiter') or '?'}  rows: {info['rows']}"
             + (f"  (+{info['bad_length_rows']} with wrong cell count)" if info.get("bad_length_rows") else "")]
    for c in info["columns"]:  # type: ignore[union-attr]
        std = f" -> {c['standard']}" if c["standard"] else ""
        ex = "; ".join(s[:40] + ("..." if len(s) > 40 else "") for s in c["samples"])
        lines.append(f"  {c['name']}{std}: {c['filled']}/{info['rows']} filled  e.g. {ex}")
    return "\n".join(lines)
