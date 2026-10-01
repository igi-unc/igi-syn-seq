#!/usr/bin/env python3
"""Recover the Kinnex MAS array profile from a real skera-segmented BAM.

Writes the adapter FASTA that `skera split` needs and a JSON profile describing the array structure, both
derived from real data rather than transcribed from a kit document. Reading them out of the data also means
the simulated arrays are segmentable by the same skera version that produced the reference set.

How the tags encode an array (measured, not assumed -- `dl` is NOT a segment count, which is how it first
looked):

    zm  the ZMW: one array
    di  the segment's index within the array, 0-7, so eight segments
    dl  the index of the adapter on the segment's LEFT
    dr  the index of the adapter on its RIGHT
    ds  MessagePack {"left": {...}, "right": {...}}, each with that adapter's `label` and observed `seq`

so the array is  A0 S0 A1 S1 A2 ... S7 A8  and a segment at index k is bracketed by adapters k and k+1.
About half of all arrays are sequenced in the reverse direction, which shows up as the mirrored pairs
(8,7), (7,6) ... (1,0) instead of (0,1), (1,2) ... (7,8).

    python3 fit_mas_adapters.py --sam <(samtools view segmented.bam) \
        --out-fasta catalog/resources/kinnex_mas8_adapters.fasta \
        --out-profile catalog/resources/kinnex_mas8_profile.json
"""
import argparse
import json
from collections import Counter, defaultdict

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_str(buf, i):
    """Read a MessagePack str at buf[i:]; (value, next) or (None, i+1)."""
    b = buf[i]
    if 0xa0 <= b <= 0xbf:
        n = b & 0x1f
        return buf[i + 1:i + 1 + n].decode("ascii", "replace"), i + 1 + n
    if b == 0xd9:
        n = buf[i + 1]
        return buf[i + 2:i + 2 + n].decode("ascii", "replace"), i + 2 + n
    return None, i + 1


def label_seq_pairs(buf):
    """Yield (label, seq) in the order the tag lists them: left adapter first, then right."""
    i, lab = 0, None
    while i < len(buf):
        s, j = read_str(buf, i)
        if s in ("label", "seq"):
            v, k = read_str(buf, j)
            if v is not None:
                if s == "label":
                    lab = v
                elif lab is not None:
                    yield lab, v
                    lab = None
                i = k
                continue
        i = i + 1 if s is None else j


def consensus(seqs):
    """Position-wise consensus over equal-length observations already in one orientation.

    Orientation is NOT resolved here. Taking each adapter's own modal sequence as its anchor gives every
    label a self-consistent but arbitrary strand, and that produced a set where adapters 0-7 were in one
    frame and adapter 8 in the other -- an array built from it would carry a reversed adapter at one end.
    The caller instead keeps only observations from forward-oriented segments, which the `dl < dr` test
    identifies, so all nine adapters arrive in the array's own forward frame.
    """
    cols = [Counter() for _ in seqs[0]]
    for s in seqs:
        for i, ch in enumerate(s):
            cols[i][ch] += 1
    seq = "".join(c.most_common(1)[0][0] for c in cols)
    worst = min(c.most_common(1)[0][1] / sum(c.values()) for c in cols)
    return seq, worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sam", required=True, help="uncompressed SAM records (samtools view output)")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-profile", required=True)
    ap.add_argument("--source", default="")
    ap.add_argument("--max-observations", type=int, default=4000,
                    help="stop decoding adapter sequences once every label has this many full-length "
                         "observations. The `ds` tag is a few hundred MessagePack bytes per segment and "
                         "decoding it in Python dominates the runtime; 4,000 observations of a 17-mer at "
                         "99.8%% accuracy settle every column far beyond doubt, so there is nothing to buy "
                         "by decoding a million of them. Structure counts (zm/di/dl/dr) are cheap and keep "
                         "being read from every record")
    a = ap.parse_args()

    obs = defaultdict(list)       # label -> observed sequences
    per_zmw = Counter()
    order = []                    # zmw ids in stream order, to drop the truncated first/last
    seglen, rq, npass = [], [], []
    grammar = Counter()
    n = 0
    enough = False

    for line in open(a.sam):
        f = line.rstrip("\n").split("\t")
        if len(f) < 11 or f[0].startswith("@"):
            continue
        ds = zm = di = dl = dr = None
        r = p = None
        for t in f[11:]:
            if t.startswith("ds:B:C,"):
                ds = bytes(int(x) & 0xff for x in t[7:].split(","))
            elif t.startswith("zm:i:"):
                zm = int(t[5:])
            elif t.startswith("di:i:"):
                di = int(t[5:])
            elif t.startswith("dl:i:"):
                dl = int(t[5:])
            elif t.startswith("dr:i:"):
                dr = int(t[5:])
            elif t.startswith("rq:f:"):
                r = float(t[5:])
            elif t.startswith("np:i:"):
                p = int(t[5:])
        if zm is None:
            continue
        n += 1
        if zm not in per_zmw:
            order.append(zm)
        per_zmw[zm] += 1
        seglen.append(len(f[9]))
        if r is not None:
            rq.append(r)
        if p is not None:
            npass.append(p)
        if None not in (di, dl, dr):
            grammar[(di, dl, dr)] += 1
        # Only forward-oriented segments (left adapter index below the right one) contribute sequence, so
        # every adapter is recorded in the array's forward frame rather than in whichever strand that read
        # happened to be sequenced on.
        if ds and not enough and dl is not None and dr is not None and dl < dr:
            for lab, v in label_seq_pairs(ds):
                obs[lab].append(v)
            # 9 adapters plus RANDOM; stop once the real ones are all well covered
            if len(obs) >= 9 and all(len(v) >= a.max_observations
                                     for k, v in obs.items() if k != "RANDOM"):
                enough = True
                print(f"  adapter sequences settled after {n:,} segments "
                      f"({min(len(v) for k, v in obs.items() if k != 'RANDOM'):,} observations of the "
                      f"least common); continuing for structure only")

    # The first and last ZMW in the stream are cut off by the head/tail of the sample, not by chemistry.
    for z in (order[0], order[-1]):
        per_zmw.pop(z, None)
    dist = Counter(per_zmw.values())
    tot_arrays = sum(dist.values())

    adapters = {}
    for lab, seqs in obs.items():
        if lab == "RANDOM":
            continue
        modal_len = Counter(len(s) for s in seqs).most_common(1)[0][0]
        full = [s for s in seqs if len(s) == modal_len]
        seq, worst = consensus(full)
        adapters[lab] = {"seq": seq, "len": modal_len, "n": len(full),
                         "min_column_agreement": round(worst, 4),
                         "frame": "array forward (from dl<dr segments only)"}

    with open(a.out_fasta, "w") as fh:
        fh.write(f"; Kinnex MAS adapters recovered from {a.source or 'a real segmented BAM'}\n"
                 f"; by jobs/measure/fit_mas_adapters.py -- see that file for how skera's tags encode them.\n")
        for lab in sorted(adapters, key=lambda x: int(x)):
            fh.write(f">{lab}\n{adapters[lab]['seq']}\n")

    sl = sorted(seglen)
    profile = {
        "source": a.source,
        "segments_per_array": 8,
        "arrays_measured": tot_arrays,
        "segments_measured": n,
        "array_yield_distribution": {str(k): round(dist[k] / tot_arrays, 5) for k in sorted(dist)},
        "mean_segments_per_array": round(sum(k * v for k, v in dist.items()) / tot_arrays, 4),
        "fraction_complete_arrays": round(dist.get(8, 0) / tot_arrays, 4),
        "segment_length": {
            "mean": round(sum(sl) / len(sl), 1),
            "median": sl[len(sl) // 2],
            "p10": sl[len(sl) // 10],
            "p90": sl[9 * len(sl) // 10],
            "sd": round((sum((x - sum(sl) / len(sl)) ** 2 for x in sl) / len(sl)) ** 0.5, 1),
        },
        "read_quality_rq_mean": round(sum(rq) / len(rq), 5) if rq else None,
        "passes_np_mean": round(sum(npass) / len(npass), 2) if npass else None,
        "fraction_arrays_reverse_oriented": round(
            sum(v for (di, dl, dr), v in grammar.items() if dl > dr) / max(1, sum(grammar.values())), 4),
        "adapters": adapters,
    }
    with open(a.out_profile, "w") as fh:
        json.dump(profile, fh, indent=2)
        fh.write("\n")

    print(f"{n:,} segments over {tot_arrays:,} arrays")
    print(f"complete (8-segment) arrays: {profile['fraction_complete_arrays']:.1%}, "
          f"mean {profile['mean_segments_per_array']:.3f} segments/array")
    print(f"segment length mean {profile['segment_length']['mean']:.0f} "
          f"sd {profile['segment_length']['sd']:.0f}, rq {profile['read_quality_rq_mean']}")
    print(f"reverse-oriented arrays: {profile['fraction_arrays_reverse_oriented']:.1%}")
    print(f"\n{len(adapters)} adapters -> {a.out_fasta}")
    for lab in sorted(adapters, key=lambda x: int(x)):
        d = adapters[lab]
        print(f"  {lab}  {d['seq']}  ({d['len']} bp, n={d['n']:,}, "
              f"min column agreement {d['min_column_agreement']:.1%})")


if __name__ == "__main__":
    main()
