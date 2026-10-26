#!/usr/bin/env python3
"""Regenerate the PDB test fixtures from the public RCSB entries.

Run once with the originals in ``~/.cache/abd-prep/`` (``abd-prep pdb fetch``
puts them there). Every fixture starts with REMARK lines that say which
entry and residues it was cut from and what was edited by hand, so the
provenance is in the file itself (plan R4). Only public coordinates are
used; nothing here is experimental or proprietary data.

    python3 tests/make_fixtures.py
"""
import os
import sys

CACHE = os.path.expanduser("~/.cache/abd-prep")
HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")


def coord_lines(pdb_id):
    with open(os.path.join(CACHE, f"{pdb_id}.pdb")) as fh:
        return [l.rstrip("\n") for l in fh if l.startswith(("ATOM", "HETATM"))]


def pick(lines, chain, lo, hi, hetatm=False):
    out = []
    for l in lines:
        if (l.startswith("HETATM") != hetatm) or l[21] != chain:
            continue
        if lo <= int(l[22:26]) <= hi:
            out.append(l)
    return out


def ter(serial, line):
    return f"TER   {serial:5d}      {line[17:20]} {line[21]}{line[22:27]}".ljust(80)


def with_altloc(line, code, occ):
    return line[:16] + code + line[17:54] + f"{occ:6.2f}" + line[60:]


def shift_x(line, dx):
    x = float(line[30:38]) + dx
    return line[:30] + f"{x:8.3f}" + line[38:]


def write(name, remarks, body):
    path = os.path.join(FIX, name)
    with open(path, "w", newline="\n") as fh:
        for i, r in enumerate(remarks, 1):
            fh.write(f"REMARK   {i} {r}".ljust(80) + "\n")
        for l in body:
            fh.write(l + "\n")
    print("wrote", path, f"({len(body)} records)")


def main():
    os.makedirs(FIX, exist_ok=True)
    n8z = coord_lines("1N8Z")
    b9z = coord_lines("6B9Z")

    # ---- 1N8Z_mini: three chains, hand-injected altlocs, waters ----------
    a = pick(n8z, "A", 1, 5)
    b = pick(n8z, "B", 1, 5)
    c = pick(n8z, "C", 489, 493)
    hoh = pick(n8z, "C", 701, 702, hetatm=True) or pick(n8z, "A", 1, 9999, hetatm=True)[:2]
    # altloc injection: CB of A:2 as A(0.60)/B(0.40); CB of B:3 as A(0.50)/B(0.50) tie
    def inject(block, resnum, occ_a, occ_b):
        out = []
        for l in block:
            if int(l[22:26]) == resnum and l[12:16].strip() == "CB":
                out.append(with_altloc(l, "A", occ_a))
                out.append(with_altloc(shift_x(l, 0.5), "B", occ_b))
            else:
                out.append(l)
        return out
    a = inject(a, 2, 0.60, 0.40)
    b = inject(b, 3, 0.50, 0.50)
    body = a + [ter(int(a[-1][6:11]) + 1, a[-1])] + b + [ter(int(b[-1][6:11]) + 1, b[-1])] \
        + c + [ter(int(c[-1][6:11]) + 1, c[-1])] + hoh + ["END".ljust(80)]
    write("1N8Z_mini.pdb", [
        "SOURCE: RCSB PDB 1N8Z (CC0) chains A 1-5, B 1-5, C 489-493 + 2 HOH",
        "EDITED: altloc A/B injected at CB of A:2 (occ 0.60/0.40) and B:3 (0.50/0.50 tie, B x+0.5)",
        "PURPOSE: altloc selection, tie rule, chain selection, TER/END regeneration",
    ], body)

    # ---- 6B9Z_cso_mini: CSO/CYS heterogeneity at B175, synthetic MSE, two MODELs
    blk = pick(b9z, "B", 173, 177) + pick(b9z, "B", 173, 177, hetatm=True)
    blk.sort(key=lambda l: int(l[6:11]))
    # synthetic MSE: turn B:177 (whatever it is) into HETATM MSE with an SE atom appended
    mse = []
    base = [l for l in blk if int(l[22:26]) == 177]
    for l in base:
        name = l[12:16].strip()
        if name in ("N", "CA", "C", "O", "CB"):
            mse.append("HETATM" + l[6:17] + "MSE" + l[20:])
    tmpl = mse[-1]
    serial = int(tmpl[6:11]) + 1
    #        1-6      7-11            12   13-16   17   18-20   21-76        77-78  79-80
    mse.append("HETATM" + f"{serial:5d}" + " " + " SE " + " " + "MSE" + tmpl[20:76] + "SE" + tmpl[78:])
    model1 = [l for l in blk if int(l[22:26]) != 177] + mse
    model2 = [shift_x(l, 1.0) for l in model1]
    body = ["MODEL        1".ljust(80)] + model1 + [ter(int(model1[-1][6:11]) + 1, model1[-1]), "ENDMDL".ljust(80),
            "MODEL        2".ljust(80)] + model2 + [ter(int(model2[-1][6:11]) + 1, model2[-1]), "ENDMDL".ljust(80),
            "END".ljust(80)]
    write("6B9Z_cso_mini.pdb", [
        "SOURCE: RCSB PDB 6B9Z (CC0) chain B 173-177 (B175 = ACSO 0.75 / BCYS 0.25 heterogeneity)",
        "EDITED: B177 rewritten as a synthetic HETATM MSE (N,CA,C,O,CB + SE); wrapped in MODEL 1 and MODEL 2 (x+1.0)",
        "PURPOSE: modified residues (CSO, MSE), --std-residues, residue heterogeneity, first-MODEL-only",
    ], body)

    # ---- kabat_mini: insertion codes 82A/82B/82C ------------------------------
    src = pick(n8z, "B", 80, 86)
    mapping = {80: (80, " "), 81: (81, " "), 82: (82, " "), 83: (82, "A"), 84: (82, "B"), 85: (82, "C"), 86: (83, " ")}
    out = []
    for l in src:
        num, code = mapping[int(l[22:26])]
        out.append(l[:22] + f"{num:4d}" + code + l[27:])
    body = out + [ter(int(out[-1][6:11]) + 1, out[-1]), "END".ljust(80)]
    write("kabat_mini.pdb", [
        "SOURCE: RCSB PDB 1N8Z (CC0) chain B 80-86",
        "EDITED: renumbered Kabat-style to 80, 81, 82, 82A, 82B, 82C, 83",
        "PURPOSE: insertion-code range selection (H:1-82 vs H:1-82A), numbering-gap rules",
    ], body)


if __name__ == "__main__":
    for pid in ("1N8Z", "6B9Z"):
        if not os.path.isfile(os.path.join(CACHE, f"{pid}.pdb")):
            sys.exit(f"missing {CACHE}/{pid}.pdb -- run: abd-prep pdb fetch {pid}")
    main()
