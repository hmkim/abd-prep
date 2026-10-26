"""Download PDB-format files and entry metadata from RCSB (F1, F9).

Two public endpoints, both HTTPS, both built from a *validated* PDB ID:

* ``https://files.rcsb.org/download/<ID>.pdb``            coordinates
* ``https://data.rcsb.org/rest/v1/core/entry/<ID>``       metadata (JSON)

Files land in a local cache (default ``~/.cache/abd-prep/``, override with
``--cache`` or ``ABD_PREP_CACHE``) and are reused on the next run. The
download is streamed to ``<ID>.pdb.part`` and renamed into place only when
complete, so an interrupted transfer never leaves a truncated file that
looks valid.

Only ``urllib`` is used; it honours ``HTTPS_PROXY`` from the environment.
``opener`` is injectable so tests run without the network.
"""
import json
import os
import re
import shutil
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

from .. import __version__
from ..common.result import InputError

RCSB_FILES_URL = "https://files.rcsb.org/download/"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/"
DEFAULT_CACHE_DIR = "~/.cache/abd-prep"
CACHE_ENV = "ABD_PREP_CACHE"
TIMEOUT_SECONDS = 60
#: Hard stop for one download. PDB-format files are at most a few MB per
#: model (99,999-atom limit); large multi-model NMR entries stay far below.
MAX_DOWNLOAD_BYTES = 512 * 1024 * 1024
_CHUNK = 256 * 1024

#: Classic 4-character PDB ID: a digit 1-9 followed by three alphanumerics.
PDB_ID_RE = re.compile(r"^[1-9][A-Za-z0-9]{3}$")

#: ABD user guide: structures of 2.5 A or better are recommended.
RECOMMENDED_RESOLUTION = 2.5

#: ``opener(url) -> (http_status, body)``. Raises :class:`FetchError` on a
#: network failure (no connection, timeout, bad certificate ...).
Opener = Callable[[str], Tuple[int, bytes]]


class FetchError(InputError):
    """Download problem; the CLI maps it to exit code 2."""


def normalize_pdb_id(text: str) -> str:
    """Validate and upper-case a PDB ID (``1n8z`` -> ``1N8Z``).

    Rejecting anything but ``[1-9][A-Za-z0-9]{3}`` is also what keeps the
    ID safe to use in a URL and as a cache file name.
    """
    text = (text or "").strip()
    if not PDB_ID_RE.match(text):
        raise InputError(f"{text!r} is not a PDB ID (expected 4 characters such as 1N8Z)")
    return text.upper()


def cache_dir(override: Optional[str] = None, create: bool = True) -> str:
    path = override or os.environ.get(CACHE_ENV) or DEFAULT_CACHE_DIR
    path = os.path.abspath(os.path.expanduser(path))
    if create:
        os.makedirs(path, exist_ok=True)
    return path


def pdb_url(pdb_id: str) -> str:
    return RCSB_FILES_URL + normalize_pdb_id(pdb_id) + ".pdb"


def entry_url(pdb_id: str) -> str:
    return RCSB_ENTRY_URL + normalize_pdb_id(pdb_id)


def default_opener(url: str) -> Tuple[int, bytes]:
    """GET ``url`` with urllib; HTTP errors become ``(status, b"")``."""
    req = urllib.request.Request(url, headers={"User-Agent": f"abd-prep/{__version__}"})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT_SECONDS) as resp:
            chunks: List[bytes] = []
            total = 0
            while True:
                chunk = resp.read(_CHUNK)
                if not chunk:
                    break
                total += len(chunk)
                if total > MAX_DOWNLOAD_BYTES:
                    raise FetchError(f"{url}: response larger than {MAX_DOWNLOAD_BYTES // (1024 * 1024)} MB; "
                                     f"refusing to continue")
                chunks.append(chunk)
            return resp.status, b"".join(chunks)
    except urllib.error.HTTPError as exc:
        return exc.code, b""
    except urllib.error.URLError as exc:
        raise FetchError(f"cannot reach {url}: {exc.reason}")
    except OSError as exc:            # socket timeout, connection reset
        raise FetchError(f"cannot reach {url}: {exc}")


@dataclass
class EntryMeta:
    """The few entry fields the report keeps (plan F9)."""

    pdb_id: str
    title: str = ""
    revision_date: str = ""
    initial_release_date: str = ""
    resolution: Optional[float] = None
    experimental_method: str = ""
    model_count: Optional[int] = None

    def to_dict(self) -> Dict[str, object]:
        return {
            "pdb_id": self.pdb_id,
            "title": self.title,
            "revision_date": self.revision_date,
            "initial_release_date": self.initial_release_date,
            "resolution_angstrom": self.resolution,
            "experimental_method": self.experimental_method,
            "model_count": self.model_count,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, object]) -> "EntryMeta":
        return cls(
            pdb_id=str(d.get("pdb_id", "")),
            title=str(d.get("title", "") or ""),
            revision_date=str(d.get("revision_date", "") or ""),
            initial_release_date=str(d.get("initial_release_date", "") or ""),
            resolution=_float_or_none(d.get("resolution_angstrom")),
            experimental_method=str(d.get("experimental_method", "") or ""),
            model_count=_int_or_none(d.get("model_count")),
        )


def summarize_entry(raw: Dict[str, object], pdb_id: str) -> EntryMeta:
    """Pick the fields we keep out of the RCSB entry JSON (structure as of 2026-10)."""
    acc = raw.get("rcsb_accession_info") or {}
    info = raw.get("rcsb_entry_info") or {}
    struct = raw.get("struct") or {}
    res_list = info.get("resolution_combined") or []
    resolution = _float_or_none(res_list[0]) if isinstance(res_list, list) and res_list else None
    methods = info.get("experimental_method") or ", ".join(
        str(e.get("method", "")) for e in (raw.get("exptl") or []) if isinstance(e, dict))
    return EntryMeta(
        pdb_id=pdb_id,
        title=str(struct.get("title", "") or ""),
        revision_date=_date_only(acc.get("revision_date")),
        initial_release_date=_date_only(acc.get("initial_release_date")),
        resolution=resolution,
        experimental_method=str(methods or ""),
        model_count=_int_or_none(info.get("deposited_model_count")),
    )


def _date_only(value: object) -> str:
    text = str(value or "")
    return text[:10] if len(text) >= 10 else text


def _float_or_none(value: object) -> Optional[float]:
    try:
        return None if value is None else float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _int_or_none(value: object) -> Optional[int]:
    try:
        return None if value is None else int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


@dataclass
class FetchResult:
    pdb_id: str
    path: str
    from_cache: bool
    bytes: int
    meta: Optional[EntryMeta] = None
    meta_error: Optional[str] = None
    notes: List[str] = field(default_factory=list)

    def to_dict(self, path_mode: str = "basename") -> Dict[str, object]:
        return {
            "pdb_id": self.pdb_id,
            "path": os.path.basename(self.path) if path_mode == "basename" else self.path,
            "from_cache": self.from_cache,
            "bytes": self.bytes,
            "entry": self.meta.to_dict() if self.meta else None,
            "entry_error": self.meta_error,
        }


def _entry_exists(pdb_id: str, opener: Opener) -> Optional[bool]:
    try:
        status, _ = opener(entry_url(pdb_id))
    except FetchError:
        return None
    return status == 200


def fetch_pdb(pdb_id: str, cache: Optional[str] = None, refresh: bool = False,
              opener: Opener = default_opener) -> FetchResult:
    """Return the cached PDB-format file for ``pdb_id``, downloading if needed."""
    pdb_id = normalize_pdb_id(pdb_id)
    directory = cache_dir(cache)
    path = os.path.join(directory, f"{pdb_id}.pdb")
    if not refresh and os.path.isfile(path) and os.path.getsize(path) > 0:
        return FetchResult(pdb_id, path, True, os.path.getsize(path))

    url = pdb_url(pdb_id)
    status, body = opener(url)
    if status == 404:
        exists = _entry_exists(pdb_id, opener)
        if exists:
            raise FetchError(f"entry {pdb_id} has no PDB-format file at RCSB (mmCIF only -- usually a structure "
                             f"too large for the legacy format). mmCIF input is not supported yet; "
                             f"see https://files.rcsb.org/download/{pdb_id}.cif")
        if exists is False:
            raise FetchError(f"entry {pdb_id} does not exist at RCSB ({url} -> 404)")
        raise FetchError(f"{url} -> 404 and the entry service could not be reached to tell why")
    if status != 200:
        raise FetchError(f"{url} -> HTTP {status}")
    if not body.lstrip().startswith((b"HEADER", b"ATOM", b"HETATM", b"REMARK", b"MODEL", b"TITLE", b"COMPND")):
        raise FetchError(f"{url}: response does not look like a PDB file")

    part = path + ".part"
    try:
        with open(part, "wb") as fh:
            fh.write(body)
        os.replace(part, path)
    finally:
        if os.path.exists(part):
            os.remove(part)
    return FetchResult(pdb_id, path, False, len(body))


def fetch_entry(pdb_id: str, cache: Optional[str] = None, refresh: bool = False,
                opener: Opener = default_opener) -> Tuple[Optional[EntryMeta], Optional[str]]:
    """Return ``(meta, error)``; metadata is best-effort and cached as ``<ID>.entry.json``."""
    pdb_id = normalize_pdb_id(pdb_id)
    directory = cache_dir(cache)
    path = os.path.join(directory, f"{pdb_id}.entry.json")
    if not refresh and os.path.isfile(path):
        try:
            with open(path, "r", encoding="utf-8") as fh:
                return EntryMeta.from_dict(json.load(fh)), None
        except (OSError, ValueError):
            pass                                  # fall through and refetch
    try:
        status, body = opener(entry_url(pdb_id))
    except FetchError as exc:
        return None, str(exc)
    if status != 200:
        return None, f"{entry_url(pdb_id)} -> HTTP {status}"
    try:
        raw = json.loads(body.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        return None, f"entry metadata is not valid JSON: {exc}"
    if not isinstance(raw, dict):
        return None, "entry metadata is not a JSON object"
    meta = summarize_entry(raw, pdb_id)
    try:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(meta.to_dict(), fh, indent=2)
            fh.write("\n")
    except OSError:
        pass                                      # cache write failures are not fatal
    return meta, None


def fetch(pdb_id: str, cache: Optional[str] = None, refresh: bool = False,
          opener: Opener = default_opener, with_meta: bool = True) -> FetchResult:
    """Coordinates plus (best-effort) entry metadata in one call."""
    result = fetch_pdb(pdb_id, cache, refresh, opener)
    if with_meta:
        result.meta, result.meta_error = fetch_entry(pdb_id, cache, refresh, opener)
    return result


def copy_to(result: FetchResult, out_dir: str) -> str:
    """Copy a fetched file into ``out_dir`` (``abd-prep pdb fetch --out``)."""
    os.makedirs(out_dir, exist_ok=True)
    dest = os.path.join(out_dir, os.path.basename(result.path))
    if os.path.abspath(dest) != os.path.abspath(result.path):
        shutil.copyfile(result.path, dest)
    return dest
