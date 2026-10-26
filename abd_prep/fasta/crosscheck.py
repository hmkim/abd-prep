"""F18: is this FASTA seed the same molecule as that PDB chain?

Use case (plan F18): the same antibody exists twice -- as an Fv structure
(e.g. the ``pdb split`` output) and as a sequence from a database or
sequencing run. Before uploading both to ABD, confirm they agree.

Why not a plain string comparison: a crystal structure has holes
(disordered residues that were not modelled), so its extracted sequence
is shorter than the real chain, and a seed FASTA may carry extra
N-/C-terminal residues (tag, signal peptide, constant region). The
comparison therefore walks the PDB chain **by residue number**: a jump in
``resSeq`` is a gap of that many residues, insertion codes (82, 82A, 82B)
count as consecutive, and the FASTA offset that best explains the chain
is searched over every possible start. Afterwards:

* PDB residues missing (gap)                 -> info     ``FASTA_PDB_GAP``
* both present but different residue         -> error    ``FASTA_PDB_MISMATCH``
* PDB residues beyond the FASTA record       -> warning  ``FASTA_PDB_UNCOVERED``
* FASTA residues beyond the PDB chain        -> info     ``FASTA_EXTRA_RESIDUES``
* best alignment still mostly different      -> error    ``FASTA_PDB_NO_CORRESPONDENCE``
* identity summary per chain                 -> info     ``FASTA_PDB_IDENTITY``
* mapped role without a FASTA record         -> warning  ``FASTA_PDB_ROLE_MISSING``
* no record could be compared at all         -> error    ``FASTA_PDB_NOT_COMPARED``

Chain mapping: when the PDB chains are already labelled ``H``/``L`` the
roles map onto them automatically; otherwise ``--chain-map H=B,L=A`` is
required (missing -> input error, exit 2).
"""
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from ..common import residues as R
from ..common.result import InputError, ValidationResult
from ..pdb.pdbfile import iter_residues, read_pdb
from .fastafile import Record
from .roles import ROLE_HEAVY, ROLE_LIGHT, ROLE_SCFV

#: Below this identity (over aligned residues) the chains are considered
#: unrelated and per-position mismatches are not listed.
MIN_IDENTITY = 0.5
#: ... unless the alignment is very short; need at least this many matches.
MIN_MATCHES = 10
#: Maximum positions listed in one mismatch finding.
LIST_LIMIT = 20


@dataclass
class ChainResidues:
    chain: str
    labels: List[str] = field(default_factory=list)          # "52", "82A" ...
    letters: str = ""                                         # one-letter sequence (gaps NOT represented)
    index: List[int] = field(default_factory=list)            # virtual sequential index, gap-aware
    gaps: List[Tuple[str, str, int, int]] = field(default_factory=list)  # (before, after, n_missing, index_after)

    @property
    def span(self) -> int:
        return (self.index[-1] + 1) if self.index else 0


def pdb_chain_residues(path: str) -> Dict[str, ChainResidues]:
    """Per-chain residue letters with gap-aware sequential indices.

    Polymer HETATM residues that ``--std-residues`` knows (MSE, CSO ...)
    are included as their parent amino acid so an un-cleaned structure does
    not show a false gap; ligands and waters are skipped.
    """
    pdb = read_pdb(path)
    out: Dict[str, ChainResidues] = {}
    prev: Dict[str, Tuple[int, str]] = {}
    for (chain, resseq, icode), resname, atoms in iter_residues(pdb.first_model_atoms()):
        if atoms[0].is_hetatm:
            mod = R.MODIFIED_RESIDUES.get(resname.upper())
            if mod is None:
                continue
            letter = R.one_letter(mod.parent)
        else:
            letter = R.one_letter(resname)
        cr = out.get(chain)
        if cr is None:
            cr = out[chain] = ChainResidues(chain)
            k = 0
        else:
            p_seq, p_icode = prev[chain]
            if resseq > p_seq:
                k = cr.index[-1] + (resseq - p_seq)
                if resseq - p_seq > 1:
                    cr.gaps.append((cr.labels[-1], f"{resseq}{icode.strip()}", resseq - p_seq - 1, k))
            else:   # same number with a new insertion code, or decreasing numbering: consecutive
                k = cr.index[-1] + 1
        cr.labels.append(f"{resseq}{icode.strip()}")
        cr.letters += letter
        cr.index.append(k)
        prev[chain] = (resseq, icode)
    return out


def parse_chain_map(text: Optional[str]) -> Dict[str, str]:
    """``"H=B,L=A"`` -> ``{"H": "B", "L": "A"}``."""
    out: Dict[str, str] = {}
    if not text:
        return out
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise InputError(f"--chain-map: expected ROLE=CHAIN (e.g. H=B), got {item!r}")
        role, chain = (p.strip() for p in item.split("=", 1))
        role = role.upper()
        if role not in (ROLE_HEAVY, ROLE_LIGHT):
            raise InputError(f"--chain-map: role must be H or L, got {role!r}")
        if len(chain) != 1:
            raise InputError(f"--chain-map: PDB chain IDs are one character, got {chain!r}")
        if role in out:
            raise InputError(f"--chain-map: role {role} given twice")
        out[role] = chain
    return out


def resolve_chain_map(chains: Dict[str, ChainResidues], chain_map: Dict[str, str], pdb_path: str) -> Dict[str, str]:
    """Role -> PDB chain. Automatic only when the PDB uses H/L labels."""
    if chain_map:
        for role, ch in chain_map.items():
            if ch not in chains:
                raise InputError(f"--chain-map {role}={ch}: the PDB has no chain {ch!r} "
                                 f"(chains: {', '.join(chains) or 'none'})")
        return dict(chain_map)
    auto = {role: role for role in (ROLE_HEAVY, ROLE_LIGHT) if role in chains}
    if not auto:
        raise InputError(f"cannot tell which PDB chain is heavy/light: {pdb_path} has chains "
                         f"{', '.join(chains) or 'none'} and none is labelled H or L. Pass --chain-map, "
                         "e.g. --chain-map H=B,L=A")
    return auto


@dataclass
class Alignment:
    offset: int = 0            # FASTA index (0-based) of PDB virtual index 0
    overlap: int = 0
    matches: int = 0
    mismatches: List[Tuple[str, str, int, str]] = field(default_factory=list)  # (label, pdb_aa, fasta_pos1, fasta_aa)
    uncovered: List[str] = field(default_factory=list)                        # PDB labels outside the record
    gaps_in_range: List[Tuple[str, str, int, int]] = field(default_factory=list)  # (before, after, n, fasta_pos1)
    extra_n: int = 0          # FASTA residues before the chain start
    extra_c: int = 0          # FASTA residues after the chain end

    @property
    def identity(self) -> float:
        return self.matches / self.overlap if self.overlap else 0.0


def align(chain: ChainResidues, seq: str) -> Alignment:
    """Best offset of the gap-aware PDB chain inside the FASTA sequence.

    Exhaustive over every offset: O(len(seq) x residues). Antibody inputs
    take milliseconds; the worst case the tool accepts (10,000-residue
    FASTA against a 1,000-residue chain) measured about 2 s.
    """
    n, idx, letters = len(seq), chain.index, chain.letters
    if not idx or not n:
        return Alignment()
    span = chain.span
    best: Optional[Tuple[int, int, int]] = None  # (matches, -mismatches, -abs(offset)) for max()
    best_o = 0
    min_overlap = min(5, len(idx))
    for o in range(-span + 1, n):
        matches = mismatches = overlap = 0
        for k, aa in zip(idx, letters):
            j = k + o
            if 0 <= j < n:
                overlap += 1
                if seq[j] == aa:
                    matches += 1
                else:
                    mismatches += 1
        if overlap < min_overlap:
            continue
        key = (matches, -mismatches, -abs(o))
        if best is None or key > best:
            best, best_o = key, o
    if best is None:
        return Alignment()
    al = Alignment(offset=best_o)
    for label, k, aa in zip(chain.labels, idx, letters):
        j = k + best_o
        if 0 <= j < n:
            al.overlap += 1
            if seq[j] == aa:
                al.matches += 1
            else:
                al.mismatches.append((label, aa, j + 1, seq[j]))
        else:
            al.uncovered.append(label)
    for before, after, n_missing, k_after in chain.gaps:
        j_first = k_after - n_missing + best_o
        if 0 <= j_first and j_first + n_missing <= n:
            al.gaps_in_range.append((before, after, n_missing, j_first + 1))
    first_j, last_j = idx[0] + best_o, idx[-1] + best_o
    al.extra_n = max(0, first_j)
    al.extra_c = max(0, n - 1 - last_j)
    return al


def _record_for(records: Sequence[Record], role: str) -> List[Record]:
    return [r for r in records if r.role == role or r.role == ROLE_SCFV]


def crosscheck(records: Sequence[Record], pdb_path: str, chain_map: Dict[str, str],
               result: ValidationResult) -> Dict[str, object]:
    """Compare every role-bearing record with its mapped PDB chain; return details."""
    chains = pdb_chain_residues(pdb_path)
    mapping = resolve_chain_map(chains, chain_map, pdb_path)
    # basename only: this entry is informational and report.json defaults to leak-free basenames;
    # the full path (when asked for with --report-paths full) is in the report's inputs list
    details: Dict[str, object] = {"pdb": os.path.basename(pdb_path), "chain_map": mapping, "chains": {}}
    compared = 0
    for role, ch in sorted(mapping.items()):
        chain = chains[ch]
        recs = _record_for(records, role)
        if not recs:
            result.warning("FASTA_PDB_ROLE_MISSING",
                           f"PDB chain {ch} is mapped to role {role} but no FASTA record has that role; not compared",
                           f"{role}:{ch}")
            continue
        for rec in recs:
            compared += 1
            where = f"{role}:{ch}<-{rec.id}" if rec.id else f"{role}:{ch}"
            al = align(chain, rec.sequence)
            details["chains"][where] = {  # type: ignore[index]
                "pdb_chain": ch, "record": rec.id, "offset": al.offset, "aligned": al.overlap,
                "identical": al.matches, "mismatches": len(al.mismatches), "uncovered_pdb": len(al.uncovered),
                "pdb_gaps": [list(g) for g in al.gaps_in_range], "extra_n": al.extra_n, "extra_c": al.extra_c,
            }
            if al.overlap == 0 or al.identity < MIN_IDENTITY or al.matches < min(MIN_MATCHES, len(chain.letters)):
                result.error("FASTA_PDB_NO_CORRESPONDENCE",
                             f"PDB chain {ch} ({len(chain.letters)} res) and record {rec.id!r} ({rec.length} res) do "
                             f"not correspond: best alignment {al.matches}/{al.overlap} identical. Wrong record, wrong "
                             "chain, or wrong --chain-map", where)
                continue
            result.info("FASTA_PDB_IDENTITY",
                        f"PDB chain {ch} vs record {rec.id!r}: {al.matches}/{al.overlap} aligned residues identical "
                        f"(PDB residue {chain.labels[0]} = FASTA position {chain.index[0] + al.offset + 1})", where)
            if al.mismatches:
                shown = "; ".join(f"PDB {lbl} {p} != FASTA {pos} {f}" for lbl, p, pos, f in al.mismatches[:LIST_LIMIT])
                more = f" ... ({len(al.mismatches)} total)" if len(al.mismatches) > LIST_LIMIT else ""
                result.error("FASTA_PDB_MISMATCH",
                             f"{len(al.mismatches)} residue(s) differ: {shown}{more}", where)
            for before, after, n_missing, pos1 in al.gaps_in_range:
                result.info("FASTA_PDB_GAP",
                            f"PDB chain {ch} has no residues between {before} and {after} ({n_missing} missing, "
                            f"FASTA positions {pos1}-{pos1 + n_missing - 1}); unmodelled in the structure", where)
            if al.uncovered:
                result.warning("FASTA_PDB_UNCOVERED",
                               f"{len(al.uncovered)} PDB residue(s) ({al.uncovered[0]}..{al.uncovered[-1]}) lie outside "
                               f"the FASTA record; the seed sequence is shorter than the structure", where)
            if al.extra_n:
                result.info("FASTA_EXTRA_RESIDUES",
                            f"FASTA has {al.extra_n} N-terminal residue(s) not in PDB chain {ch} "
                            "(signal peptide or tag?)", where)
            if al.extra_c:
                result.info("FASTA_EXTRA_RESIDUES",
                            f"FASTA has {al.extra_c} C-terminal residue(s) not in PDB chain {ch} "
                            "(constant region or tag?)", where)
    if compared == 0:
        result.error("FASTA_PDB_NOT_COMPARED",
                     "--against compared nothing: no FASTA record carries a heavy/light role. Give --framework "
                     "(nanobody/scfv) or --heavy-id/--light-id so the records can be matched to the PDB chains",
                     "against")
    details["compared"] = compared
    return details
