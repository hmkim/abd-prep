"""Clean a PDB structure for ABD upload (F3, F4, F5, F11, F12).

This is the prototype ``clean_pdb.py`` turned into a library:

1. **altloc pass** -- for every atom that has alternate locations, keep
   the one with the highest occupancy (first seen on a tie) and blank
   column 17. Same algorithm and tie rule as the prototype, so outputs
   are byte-identical for the four reference files.
2. **filter pass** -- first MODEL only (P4), HETATM removal with
   modified-residue protection (P3), chain/residue selection with
   insertion codes (P2), optional hydrogen stripping, chain renaming.

TER/END regeneration (P1) happens in :func:`abd_prep.pdb.pdbfile.write_pdb`.
Nothing here changes coordinates, B-factors or atom order.
"""
from collections import Counter, OrderedDict
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Set, Tuple

from ..common import residues as R
from ..common.result import ValidationResult
from .pdbfile import Atom, PdbFile, iter_residues, read_pdb, ter_serial_collisions, write_pdb
from .ranges import ChainSelection

AltKey = Tuple[str, str, str, str]           # (chain, resSeq text, iCode, atom-name field)
ResKey = Tuple[str, int, str]                # (chain, resSeq, iCode)
BACKBONE_TRIAD = frozenset({"N", "CA", "C"})


@dataclass
class CleanOptions:
    selection: ChainSelection = field(default_factory=ChainSelection)
    rename: Dict[str, str] = field(default_factory=dict)
    keep_het: bool = False          # prototype --keep-het
    std_residues: bool = False      # F4: convert CSO/MSE/... to the standard parent
    strip_h: bool = False           # F12
    renumber: bool = False          # rewrite atom/TER serials 1..N on output


@dataclass
class CleanOutcome:
    atoms: List[Atom]
    result: ValidationResult
    details: Dict[str, object]
    ters: List[str] = field(default_factory=list)


def select_altlocs(atoms: List[Atom], candidate: Optional[Callable[[Atom], bool]] = None) -> Dict[AltKey, Tuple[float, str]]:
    """Pass 1 of the prototype: best (occupancy, altloc) per atom key.

    ``candidate`` restricts the contest to atoms that will actually be
    written. Without it a HETATM conformer (e.g. CSO altloc A, occ 0.75)
    could win the contest and then be removed as HETATM, taking the
    surviving ATOM conformer (CYS altloc B, occ 0.25) down with it -- the
    residue would vanish. Seen at 6B9Z B175.
    """
    best: Dict[AltKey, Tuple[float, str]] = {}
    for a in atoms:
        alt = a.altloc
        if alt == " ":
            continue
        if candidate is not None and not candidate(a):
            continue
        key = (a.line[21], a.line[22:26], a.line[26], a.line[12:16])
        occ = a.occupancy
        occ = 0.0 if occ is None else occ
        if key not in best or occ > best[key][0]:
            best[key] = (occ, alt)
    return best


def polymer_like_hetero_residues(atoms: List[Atom]) -> Dict[ResKey, str]:
    """HETATM residues that look like part of the protein chain (P3).

    A HETATM residue is "polymer-like" when it carries the backbone atoms
    N, CA and C and its chain also has ATOM records. Water, glycans, ions
    and most ligands fail the backbone test; CSO, MSE, SEP, PCA ... pass.
    Returns ``{residue key: resname}``.
    """
    chains_with_atom: Set[str] = {a.chain for a in atoms if not a.is_hetatm}
    names: Dict[ResKey, Set[str]] = {}
    resnames: Dict[ResKey, str] = {}
    for a in atoms:
        if not a.is_hetatm:
            continue
        key = a.residue_key
        names.setdefault(key, set()).add(a.name)
        resnames[key] = a.resname
    return {k: resnames[k] for k, ns in names.items()
            if BACKBONE_TRIAD <= ns and k[0] in chains_with_atom}


def residue_winners(atoms: List[Atom], candidate: Optional[Callable[[Atom], bool]] = None) -> Dict[ResKey, Tuple[str, Dict[str, float]]]:
    """Resolve residue-level heterogeneity (two residue *names* at one position).

    Some entries model e.g. ``ACSO``/``BCYS`` or ``ASER``/``BGLY`` at the same
    (chain, resSeq, iCode). Picking altlocs atom by atom would then mix
    atoms of two residue types into one corrupt residue, so the residue
    name with the highest occupancy wins first and every atom of the other
    names is dropped. Returns ``{residue key: (winner, {name: max occ})}``
    for the affected residues only.
    """
    occ_by_name: Dict[ResKey, Dict[str, float]] = {}
    for a in atoms:
        if a.altloc == " " or (candidate is not None and not candidate(a)):
            continue
        occ = a.occupancy
        occ = 0.0 if occ is None else occ
        names = occ_by_name.setdefault(a.residue_key, {})
        names[a.resname] = max(occ, names.get(a.resname, -1.0))
    out: Dict[ResKey, Tuple[str, Dict[str, float]]] = {}
    for key, names in occ_by_name.items():
        if len(names) > 1:
            winner = max(names.items(), key=lambda kv: kv[1])[0]   # first-seen wins a tie (dict order)
            out[key] = (winner, names)
    return out


def clean_atoms(pdb: PdbFile, opts: CleanOptions) -> CleanOutcome:
    """Run both passes over a parsed file and return the output atoms."""
    result = ValidationResult()
    atoms = pdb.first_model_atoms()
    if pdb.n_models > 1:
        result.warning("PDB_MULTI_MODEL",
                       f"{pdb.n_models} MODEL blocks; only the first model is kept", "file")

    polymer_like = polymer_like_hetero_residues(atoms)
    sel = opts.selection

    def conversion_rule(a: Atom) -> Optional[R.ModifiedResidue]:
        if not (opts.std_residues and a.is_hetatm):
            return None
        rule = R.MODIFIED_RESIDUES.get(a.resname)
        return rule if rule is not None and a.residue_key in polymer_like else None

    def survives_het_filter(a: Atom) -> bool:
        return (not a.is_hetatm) or opts.keep_het or conversion_rule(a) is not None

    best = select_altlocs(atoms, survives_het_filter)
    winners = residue_winners(atoms, survives_het_filter)
    atom_residue_keys: Set[ResKey] = {a.residue_key for a in atoms if not a.is_hetatm}

    # -- selection sanity: every requested chain must exist ---------------
    present_chains: Set[str] = {a.chain for a in atoms}
    for ch in sel.chains:
        if ch not in present_chains:
            result.error("PDB_CHAIN_NOT_FOUND",
                         f"chain {ch!r} is not in the file (chains present: {''.join(sorted(present_chains))})",
                         ch)

    out: List[Atom] = []
    altloc_dropped = 0
    altloc_resolved_atoms = 0
    altloc_resolved_res: Set[ResKey] = set()
    het_removed: Counter = Counter()
    hydrogens_removed = 0
    modified_warned: Set[ResKey] = set()
    converted: "OrderedDict[ResKey, Dict[str, object]]" = OrderedDict()
    conv_dropped: Dict[ResKey, List[str]] = {}
    out_chain_sources: Dict[str, Set[str]] = {}
    selected_sources: Set[str] = set()
    het_only_sources: Set[str] = set()        # chains whose in-range atoms were all removed HETATM
    non_ascii = 0

    for a in atoms:
        chain, resseq, icode, resname = a.chain, a.resseq, a.icode, a.resname
        rkey = (chain, resseq, icode)
        in_range = sel.wants(chain, resseq, icode)
        where = f"{chain}:{resseq}{icode.strip()}"

        # 1. HETATM filter (F4 / P3) -------------------------------------
        if a.is_hetatm and not survives_het_filter(a):
            het_removed[resname] += 1
            if in_range:
                het_only_sources.add(chain)
            if in_range and rkey in polymer_like and rkey not in modified_warned:
                modified_warned.add(rkey)
                rule = R.MODIFIED_RESIDUES.get(resname)
                if rkey in atom_residue_keys:
                    result.info("PDB_MODIFIED_CONFORMER_REMOVED",
                                f"{resname} conformer removed with the HETATM records; the standard-residue "
                                f"conformer at this position is kept, no gap", where)
                elif rule is not None:
                    result.warning("PDB_MODIFIED_RESIDUE_REMOVED",
                                   f"{resname} is a modified backbone residue inside the selection and was "
                                   f"removed with the other HETATM records -> gap in the chain; "
                                   f"use --std-residues to convert it to {rule.parent}", where)
                else:
                    result.warning("PDB_MODIFIED_RESIDUE_NO_RULE",
                                   f"{resname} is a modified backbone residue inside the selection and was "
                                   f"removed; no conversion rule exists (use --keep-het to retain it)", where)
            continue

        # 2. chain / residue range (F2 / P2) -------------------------------
        if not in_range:
            continue
        atom = a.copy()

        # 3. altloc: residue-level winner first, then per-atom occupancy (F3)
        alt = atom.altloc
        if alt != " ":
            if rkey in winners and resname != winners[rkey][0]:
                altloc_dropped += 1
                continue
            key = (atom.line[21], atom.line[22:26], atom.line[26], atom.line[12:16])
            if best.get(key, (0.0, alt))[1] != alt:
                altloc_dropped += 1
                continue
            atom.set_altloc(" ")
            altloc_resolved_atoms += 1
            altloc_resolved_res.add(rkey)

        # 4. modified residue -> standard parent (F4, --std-residues) --------
        rule = conversion_rule(a)
        if rule is not None:
            new_name = rule.rename.get(atom.name, atom.name)
            converted.setdefault(rkey, {"from": resname, "to": rule.parent, "note": rule.note})
            if new_name not in R.STANDARD_ATOMS[rule.parent]:
                conv_dropped.setdefault(rkey, []).append(atom.name)
                continue
            atom.set_record("ATOM")
            atom.set_resname(rule.parent)
            if new_name != atom.name:
                atom.set_name(new_name)
                el = R.element_for_atom(new_name)
                if el:
                    atom.set_element(el)

        # 5. hydrogens (F12) ---------------------------------------------------
        if opts.strip_h and atom.is_hydrogen:
            hydrogens_removed += 1
            continue

        # 6. chain rename (F5), written last ---------------------------------
        new_chain = opts.rename.get(chain, chain)
        if new_chain != chain:
            atom.set_chain(new_chain)
        if not atom.is_ascii:
            non_ascii += 1
        out_chain_sources.setdefault(new_chain, set()).add(chain)
        selected_sources.add(chain)
        out.append(atom)

    # -- findings ----------------------------------------------------------
    for rkey, (winner, names) in winners.items():
        if sel.wants(*rkey):
            listing = "/".join(f"{n}({occ:.2f})" for n, occ in names.items())
            result.info("PDB_RESIDUE_HETEROGENEITY",
                        f"two residue types modelled at one position: {listing}; kept {winner}",
                        f"{rkey[0]}:{rkey[1]}{rkey[2].strip()}")
    if altloc_resolved_atoms:
        result.info("PDB_ALTLOC_RESOLVED",
                    f"{altloc_resolved_atoms} atoms in {len(altloc_resolved_res)} residues had alternate "
                    f"locations; kept the highest-occupancy conformer, dropped {altloc_dropped} atoms", "file")
    if het_removed:
        summary = ", ".join(f"{name} x{n}" for name, n in sorted(het_removed.items()))
        result.info("PDB_HETATM_REMOVED", f"removed HETATM records: {summary}", "file")
    for rkey, info in converted.items():
        dropped = conv_dropped.get(rkey, [])
        extra = f"; dropped atoms {','.join(dropped)}" if dropped else ""
        result.info("PDB_MODIFIED_RESIDUE_CONVERTED",
                    f"{info['from']} -> {info['to']}: {info['note']}{extra}",
                    f"{rkey[0]}:{rkey[1]}{rkey[2].strip()}")
    if hydrogens_removed:
        result.info("PDB_HYDROGENS_REMOVED", f"removed {hydrogens_removed} hydrogen atoms", "file")
    for ch in sel.chains:
        if ch in present_chains and ch not in selected_sources:
            if ch in het_only_sources:
                result.error("PDB_SELECTION_EMPTY",
                             f"chain {ch!r}: the requested range (ranges: {sel.describe().get(ch)}) matched only "
                             f"HETATM records, which were removed; use --keep-het to retain them", ch)
            else:
                result.error("PDB_SELECTION_EMPTY",
                             f"chain {ch!r} is present but the requested range selected no atoms "
                             f"(ranges: {sel.describe().get(ch)})", ch)
    for new_chain, sources in out_chain_sources.items():
        if len(sources) > 1:
            result.error("PDB_RENAME_COLLISION",
                         f"chains {''.join(sorted(sources))} all map to output chain {new_chain!r}; "
                         f"they would merge into one chain", new_chain)
    if non_ascii:
        result.warning("PDB_NON_ASCII",
                       f"{non_ascii} ATOM/HETATM record(s) contain non-ASCII bytes; they are written unchanged "
                       f"but PDB is an ASCII format and parsers may reject them", "file")
    if not out:
        result.error("PDB_NO_ATOMS", "no atoms left after cleaning", "file")

    details = _details(pdb, atoms, out, altloc_resolved_atoms, altloc_dropped, het_removed,
                       hydrogens_removed, converted, opts)
    return CleanOutcome(out, result, details)


def _details(pdb: PdbFile, atoms_in: List[Atom], out: List[Atom], resolved: int, dropped: int,
             het_removed: Counter, hydrogens: int, converted, opts: CleanOptions) -> Dict[str, object]:
    per_chain: "OrderedDict[str, Dict[str, object]]" = OrderedDict()
    for key, resname, res_atoms in iter_residues(out):
        ch = key[0]
        d = per_chain.setdefault(ch, {"residues": 0, "atoms": 0, "first": None, "last": None})
        d["residues"] += 1
        d["atoms"] += len(res_atoms)
        label = f"{key[1]}{key[2].strip()}"
        if d["first"] is None:
            d["first"] = label
        d["last"] = label
    return {
        "models_in_file": pdb.n_models,
        "atoms_in": len(atoms_in),
        "atoms_out": len(out),
        "altloc_atoms_resolved": resolved,
        "altloc_atoms_dropped": dropped,
        "hetatm_removed": dict(sorted(het_removed.items())),
        "hydrogens_removed": hydrogens,
        "modified_residues_converted": [
            {"where": f"{k[0]}:{k[1]}{k[2].strip()}", "from": v["from"], "to": v["to"]} for k, v in converted.items()
        ],
        "selection": opts.selection.describe(),
        "rename": dict(opts.rename),
        "chains_out": per_chain,
    }


def clean_file(in_path: str, out_path: str, opts: CleanOptions) -> CleanOutcome:
    """Read ``in_path``, clean, write ``out_path`` with fresh TER/END."""
    pdb = read_pdb(in_path)
    outcome = clean_atoms(pdb, opts)
    outcome.ters = write_pdb(outcome.atoms, out_path, renumber=opts.renumber)
    outcome.details["ter_records_written"] = len(outcome.ters)
    collisions = [] if opts.renumber else ter_serial_collisions(outcome.atoms, outcome.ters)
    if collisions:
        shown = ", ".join(str(s) for s in collisions[:5]) + (" ..." if len(collisions) > 5 else "")
        outcome.result.warning("PDB_SERIAL_DUPLICATE",
                               f"regenerated TER serial(s) {shown} duplicate an atom serial because the source "
                               f"had no TER at that chain boundary; use --renumber for unique serials", "file")
    return outcome
