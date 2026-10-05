#!/usr/bin/env python3
"""Build the tracked germline-junction resource: per V the CDR3 nucleotides from the conserved Cys, per J
the CDR3 nucleotides up to the conserved Phe, plus the framework either side."""
import gzip, json, os, re, sys, collections
sys.path.insert(0, ".")
import yaml
from igi_catalog.genome import Genome

COMP = str.maketrans("ACGTN", "TGCAN")
_NAME = re.compile(r'gene_name "(TR[AB][VDJC][^"]*)"')
_TYPE = re.compile(r'gene_type "([^"]+)"')
B = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TAB = {a + b + c: AAS[i * 16 + j * 4 + k]
       for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}

def tr(s):
    return "".join(TAB.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))

paths = yaml.safe_load(open(sys.argv[1]))
g = Genome(paths["reference_fasta"])
exons, gtype = collections.defaultdict(set), {}
with gzip.open(paths["gtf"], "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = _NAME.search(f[8])
        if not m:
            continue
        exons[m.group(1)].add((f[0], int(f[3]) - 1, int(f[4]), f[6]))
        t = _TYPE.search(f[8])
        if t:
            gtype[m.group(1)] = t.group(1)

segs = {}
for name, ex in exons.items():
    ex = sorted(ex)
    s = "".join(g.seq(c, a, b) for c, a, b, _ in ex)
    if ex[0][3] == "-":
        s = s.translate(COMP)[::-1]
    if len(s) >= 30:
        segs[name] = s.upper()

out = {"_provenance": f"GENCODE {os.path.basename(paths['gtf'])} TR segments over GRCh38, exon sequence "
                      "joined and reverse-complemented as the strand requires",
       "_why": ("CDR3s were random peptides: `CAS` or `CAV`, then uniform draws over all twenty amino "
                "acids, then `F`. That put an internal cysteine in about half of them "
                "(CASCAWCQESIIYPATQVF) and, more importantly, made the junction independent of the V and "
                "J genes the truth table named, so a caller recovering the real V gene would disagree "
                "with the truth. Junctions are now built from these germline anchors."),
       "v": {}, "j": {}, "c": {}}

nv = nj = 0
for name in sorted(segs):
    s, ty = segs[name], gtype.get(name, "")
    if re.match(r"TR[AB]V", name) and ty == "TR_V_gene":
        best = None
        for f in range(3):
            p = tr(s[f:])
            for m in re.finditer(r"[YFH][YFLHVIAC]C", p):
                cys = f + (m.end() - 1) * 3
                if cys >= len(s) - 40 and "*" not in p[:m.end()]:
                    best = cys
        if best is None:
            continue
        # keep the CDR3 nucleotides in frame: drop a dangling 1-2 nt at the 3' end
        tail = s[best:]
        tail = tail[:len(tail) - len(tail) % 3]
        if len(tail) < 3 or not tr(tail).startswith("C"):
            continue
        out["v"][name] = {"framework": s[:best], "cdr3_nt": tail}
        nv += 1
    elif re.match(r"TR[AB]J", name) and ty == "TR_J_gene":
        got = None
        for f in range(3):
            p = tr(s[f:])
            m = re.search(r"FG[A-Z]G", p)
            if m and "*" not in p[:m.start()]:
                got = (f, f + m.start() * 3)
        if got is None:
            continue
        f, phe = got
        out["j"][name] = {"cdr3_nt": s[f:phe + 3], "framework": s[phe + 3:]}
        nj += 1
    elif re.match(r"TR[AB]C", name) and ty == "TR_C_gene":
        out["c"][name] = {"framework": s}

p = "resources/vdj_germline.json"
with open(p, "w") as fh:
    json.dump(out, fh, indent=1)
print(f"  V genes {nv}  J genes {nj}  C genes {len(out['c'])}  -> {p} ({os.path.getsize(p)} bytes)")
for fam in ("TRBV", "TRBJ", "TRAV", "TRAJ"):
    k = "v" if fam[3] == "V" else "j"
    ns = [n for n in out[k] if n.startswith(fam)]
    print(f"  {fam}: {len(ns)}  e.g. {ns[0]} cdr3_nt={out[k][ns[0]]['cdr3_nt']} -> {tr(out[k][ns[0]]['cdr3_nt'])}")
