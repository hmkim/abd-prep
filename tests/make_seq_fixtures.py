"""Generate the FASTA / CSV test fixtures (M2.5) as exact bytes.

Run from the repository root::

    python3 tests/make_seq_fixtures.py          # writes tests/fixtures/seed_*.fasta, assay_*.csv
    python3 tests/make_seq_fixtures.py --check  # exit 1 if any committed fixture differs

The fixtures are committed; this script exists so their provenance is
reproducible and so byte-level features (CR/LF line ends, a UTF-8 BOM,
CP949 encoding) are produced deterministically instead of by hand.

Provenance (plan R4: public sequences, synthetic numbers only):

* VL / VH = trastuzumab Fv from RCSB entry 1N8Z (CC0), chain A residues
  1-107 and chain B residues 1-120, extracted with
  ``abd-prep pdb inspect --fasta`` from the ``pdb clean`` output.
* The "nanobody" is SYNTHETIC: the 1N8Z VH with the four camelid
  hallmark substitutions V37F / G44E / L45R / W47G (this VH has no
  insertion before position 52, so sequential = Kabat numbering here).
  It is a test pattern for the VHH role/length rules, not a real VHH.
* The scFv is VH-(G4S)3-VL, a synthetic construct.
* The leader peptide and the C-terminal extension in the full-length
  fixture are synthetic strings with the hydrophobic / non-repetitive
  properties the heuristics look for; they are not natural sequences.
* Every number in the CSV fixtures is made up.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(HERE, "fixtures")

VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSG"
      "TDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
      "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
assert len(VL) == 107 and len(VH) == 120

SRC_VL = "RCSB 1N8Z chain A residues 1-107 (trastuzumab VL; public, CC0)"
SRC_VH = "RCSB 1N8Z chain B residues 1-120 (trastuzumab VH; public, CC0)"


def mutate(seq, *subs):
    """Apply ``(position_1based, from, to)`` substitutions, checking the original letter."""
    s = list(seq)
    for pos, old, new in subs:
        assert s[pos - 1] == old, (pos, s[pos - 1], old)
        s[pos - 1] = new
    return "".join(s)


# synthetic VHH: camelid hallmarks at 37/44/45/47
VHH = mutate(VH, (37, "V", "F"), (44, "G", "E"), (45, "L", "R"), (47, "W", "G"))
LINKER = "GGGGS" * 3
SCFV = VH + LINKER + VL
LEADER = "MGWSLILLFLVAVATRVLS"                      # synthetic: M + hydrophobic core (19 aa)
CTERM = ("ASPQTKEPLADNSVKTEYQRPLSEAVDTNGKQTVPSEDRKLTAESNPVQGDTKRSYEAPLDKNSTVQE")  # synthetic, 68 aa
assert len(CTERM) == 68 and "GGGS" not in CTERM and "GGGGS" not in CTERM
DNA = "ATGGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCTAGCATCGATCGATCGATGCTAGCTAGCATCGATCG"  # synthetic pattern


def chunks(seq, width):
    return [seq[i:i + width] for i in range(0, len(seq), width)]


def fasta_text(records, eol="\n", width=0):
    out = []
    for header, seq in records:
        out.append(">" + header)
        out.extend(chunks(seq, width) if width else [seq])
    return eol.join(out) + eol


def build():
    files = {}

    # seed_igg.fasta: header variants, CRLF, lower-case, gaps, trailing '*'
    vh_messy = VH[:66] + "\r\n" + VH[66:80].lower() + "--" + VH[80:]          # gap run inside, lower-case run
    vl_messy = VL[:60] + "\r\n" + VL[60:] + "*"
    files["seed_igg.fasta"] = (
        ">heavy trastuzumab VH, " + SRC_VH + " -- fixture: CRLF, lower-case, gaps\r\n" + vh_messy + "\r\n"
        ">VL_1 kappa trastuzumab VL, " + SRC_VL + " -- fixture: trailing stop\r\n" + vl_messy + "\r\n"
    ).encode("ascii")

    # seed_nanobody.fasta: single VHH record (synthetic, see module docstring)
    files["seed_nanobody.fasta"] = fasta_text([
        ("VHH_synthetic nanobody-like VHH (SYNTHETIC: " + SRC_VH + " with camelid hallmark substitutions "
         "V37F/G44E/L45R/W47G; not a natural sequence)", VHH),
    ], width=60).encode("ascii")

    # seed_scfv.fasta: VH-(G4S)3-VL in one record
    files["seed_scfv.fasta"] = fasta_text([
        ("scFv_1 VH-(G4S)3-VL (synthetic linker; VH " + SRC_VH + "; VL " + SRC_VL + ")", SCFV),
    ], width=60).encode("ascii")

    # seed_bad_alphabet.fasta: X/B in heavy + internal '*', Z in light
    vh_bad = mutate(VH, (10, "G", "X"), (50, "R", "B"))
    vh_bad = vh_bad[:60] + "*" + vh_bad[60:]
    vl_bad = mutate(VL, (20, "T", "Z"))
    files["seed_bad_alphabet.fasta"] = fasta_text([
        ("heavy_bad " + SRC_VH + " -- fixture: X@10, B@50, internal '*' after 60", vh_bad),
        ("light_bad " + SRC_VL + " -- fixture: Z@20", vl_bad),
    ]).encode("ascii")

    # seed_bad_nucleic.fasta: a DNA-looking heavy record + a normal light record
    files["seed_bad_nucleic.fasta"] = fasta_text([
        ("heavy_dna synthetic nucleotide pattern -- fixture: looks like DNA", DNA),
        ("light " + SRC_VL, VL),
    ]).encode("ascii")

    # seed_bad_count.fasta: three records for --framework igg (two heavy)
    files["seed_bad_count.fasta"] = fasta_text([
        ("heavy " + SRC_VH, VH),
        ("light " + SRC_VL, VL),
        ("VHH_extra synthetic VHH (see seed_nanobody.fasta) -- fixture: third record", VHH),
    ]).encode("ascii")

    # seed_warn_fulllength.fasta: leader + VH + synthetic C-terminal extension; light chain plain
    files["seed_warn_fulllength.fasta"] = fasta_text([
        ("heavy_full SYNTHETIC leader (" + str(len(LEADER)) + " aa) + " + SRC_VH + " + synthetic C-terminal "
         "extension (" + str(len(CTERM)) + " aa) -- fixture: signal-peptide and length warnings",
         LEADER + VH + CTERM),
        ("light_full " + SRC_VL, VL),
    ], width=60).encode("ascii")

    # assay_ok.csv: standard header, 20 synthetic variants, no warnings
    rows = ["id,sequence,value,unit,assay"]
    for i in range(1, 21):
        pos = 30 + i                                   # positions 31..50, inside CDR-H1/FR2/CDR-H2 region
        var = VH[:pos - 1] + ("A" if VH[pos - 1] != "A" else "G") + VH[pos:]
        kd = round(0.4 + 0.37 * i, 2)                  # synthetic KD values 0.77 .. 7.80
        rows.append(f"Ab{i:03d},{var},{kd},nM,SPR")
    files["assay_ok.csv"] = ("\n".join(rows) + "\n").encode("utf-8")

    # assay_messy_utf8bom.csv: BOM, ';', missing tokens, censored, decimal comma, mixed units, duplicate id
    rows = [
        "Clone;AA_seq;KD;Unit;Note",
        f"Ab001;{VH};2.3;nM;baseline",
        f"Ab002;{VH[:40] + 'A' + VH[41:]};NA;nM;no measurement",
        f"Ab003;{VH[:50] + 'G' + VH[51:]};<0.1;nM;below assay limit",
        f"Ab004;{VH[:55] + 'S' + VH[56:]};1,2;nM;decimal comma",
        f"Ab005;{VH[:60] + 'T' + VH[61:]};850;pM;unit differs",
        f"Ab001;{VH};2.4;nM;duplicate id",
        f"Ab006;;5.1;nM;missing sequence",
        f"Ab007;{VH[:70].lower() + VH[70:]};;nM;missing value, lower-case sequence",
    ]
    files["assay_messy_utf8bom.csv"] = b"\xef\xbb\xbf" + ("\r\n".join(rows) + "\r\n").encode("utf-8")

    # assay_messy_cp949.csv: CP949 (no BOM), Korean header, one row with an extra cell
    rows = [
        "클론,서열,KD,단위,비고",
        f"Ab001,{VH},2.3,nM,기준 서열",
        f"Ab002,{VH[:40] + 'A' + VH[41:]},3.1,nM,변이체 1",
        f"Ab003,{VH[:50] + 'G' + VH[51:]},0.9,nM,변이체 2, 추가 메모",   # unquoted comma -> 6 cells
        f"Ab004,{VH[:55] + 'S' + VH[56:]},4.4,nM,변이체 3",
    ]
    files["assay_messy_cp949.csv"] = ("\n".join(rows) + "\n").encode("cp949")
    return files


def main(argv):
    check = "--check" in argv
    files = build()
    os.makedirs(FIXTURES, exist_ok=True)
    differ = []
    for name, data in sorted(files.items()):
        path = os.path.join(FIXTURES, name)
        if check:
            with open(path, "rb") as fh:
                if fh.read() != data:
                    differ.append(name)
            continue
        with open(path, "wb") as fh:
            fh.write(data)
        print(f"wrote {path} ({len(data)} bytes)")
    if check:
        if differ:
            print("fixtures differ from the generator: " + ", ".join(differ))
            return 1
        print(f"{len(files)} fixtures match the generator")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
