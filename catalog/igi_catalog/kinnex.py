"""Kinnex (MAS-seq) array construction: full-length cDNA molecules concatenated for PacBio sequencing.

A Kinnex library ligates several full-length cDNAs into one long molecule separated by known adapters, and
sequences that as a single HiFi read. `skera split` then cuts the read back into segments. So the simulation
has to build the array, not the segments: it concatenates molecules with the real adapters, hands each array
to the read simulator as one insert, and lets the real skera do the splitting. Segmenting it ourselves would
test nothing -- the point of shipping the pre-skera BAM is that skera's own output is reproduced.

The array grammar, measured from a real segmented BAM by jobs/measure/fit_mas_adapters.py:

    array = A0 S0 A1 S1 A2 S2 A3 S3 A4 S4 A5 S5 A6 S6 A7 S7 A8

Eight segments, nine adapters, and a segment at index k is bracketed by adapters k and k+1. Adapters 0-7 are
17 bp and adapter 8 is 14 bp, all recovered in the array's forward frame by taking them only from
forward-oriented segments. About 53 % of real arrays are sequenced in the reverse direction, which shows up
in skera's tags as mirrored adapter pairs (8,7),(7,6)... instead of (0,1),(1,2)...; that is supplied by the
read simulator choosing a strand per read, so arrays are built in the forward frame only.

Two modelling decisions that matter more than they look:

1. **Molecules are sampled by abundance alone, not abundance x length.** Kinnex yields one read per
   molecule, and TPM is already length-normalised. The short-read weight is abundance x length, and on the
   release manifest the two give mean cDNA lengths of 1,788 bp and 4,197 bp -- a 2.3x error in which
   transcripts dominate the library.
2. **Size selection is applied.** Even with the right weight, simulated molecules are shorter than real
   segments (p10 498 bp against 1,355). Library prep loses short cDNA, and that loss is fitted as an
   efficiency per length bin by jobs/measure/fit_kinnex_size_selection.py.

An array carrying fewer than eight segments is built with fewer molecules, drawn from the measured yield
distribution: 77.4 % of real arrays give the full eight. That models incompleteness as failed ligation
rather than as a truncated read, which is a choice -- skera failing to find an adapter would look the same
in the output -- and it is the one that keeps the simulated segment count controllable.
"""
import json
import os
import random

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    """Minimal FASTA reader; skips the `;` comment lines our adapter file carries."""
    out, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


class SizeSelection:
    """cDNA length -> relative chance of reaching a read, from a fitted curve.

    Disabled (every length equally likely) when no curve is configured, so the builder runs before the
    curve has been fitted rather than refusing to start.
    """

    def __init__(self, path=None):
        self.bins = []          # (lo, hi, efficiency)
        self.enabled = False
        if not path or not os.path.exists(path):
            return
        d = json.load(open(path))
        for lab, eff in d["efficiency_by_length_bin"].items():
            lo, hi = lab.split("-")
            self.bins.append((int(lo), float("inf") if hi == "inf" else int(hi), float(eff)))
        self.bins.sort()
        self.enabled = bool(self.bins)

    def factor(self, length):
        if not self.enabled:
            return 1.0
        for lo, hi, e in self.bins:
            if lo <= length < hi:
                return e
        # outside every fitted bin: the curve says nothing, so do not invent an efficiency
        return 0.0


class MasArrays:
    """Build Kinnex arrays from transcript records."""

    def __init__(self, adapters_fasta, profile_json, size_selection=None):
        self.adapters = read_fasta(adapters_fasta)
        self.profile = json.load(open(profile_json))
        self.n_seg = int(self.profile.get("segments_per_array", 8))
        if len(self.adapters) != self.n_seg + 1:
            raise ValueError(f"{len(self.adapters)} adapters for {self.n_seg} segments; "
                             f"an n-segment array needs n+1 adapters")
        self.order = sorted(self.adapters, key=int)
        yd = self.profile.get("array_yield_distribution") or {}
        self.yield_k = [int(k) for k in sorted(yd, key=int)]
        self.yield_w = [float(yd[str(k)]) for k in self.yield_k]
        # Recorded, not applied: the read simulator randomises strand per read, so this is a target to
        # check the built library against rather than a knob to turn.
        self.p_reverse_expected = float(self.profile.get("fraction_arrays_reverse_oriented", 0.5))
        self.size = size_selection or SizeSelection(None)

    def sample_molecules(self, records, n, rng):
        """Draw n molecules with replacement, weighted by abundance and size-selection efficiency."""
        live = [r for r in records if r.get("abundance", 0) > 0 and r.get("sequence")]
        if not live:
            raise ValueError("no records with positive abundance and a sequence")
        w = [r["abundance"] * self.size.factor(len(r["sequence"])) for r in live]
        if sum(w) <= 0:
            raise ValueError("every molecule weight is zero; the size-selection curve excludes everything")
        return rng.choices(live, weights=w, k=n)

    def segments_per_array(self, rng):
        if not self.yield_k:
            return self.n_seg
        return rng.choices(self.yield_k, weights=self.yield_w, k=1)[0]

    def build(self, molecules, rng, start_index=0):
        """Group molecules into arrays. Yields (array_name, sequence, [member records]).

        Molecules are consumed in order, so the caller controls their composition; the array size is drawn
        per array from the measured yield distribution.
        """
        i = 0
        idx = start_index
        while i < len(molecules):
            k = min(self.segments_per_array(rng), self.n_seg, len(molecules) - i)
            if k < 1:
                break
            members = molecules[i:i + k]
            i += k
            parts = []
            for j, m in enumerate(members):
                parts.append(self.adapters[self.order[j]])
                parts.append(self._random_spacer(rng))
                parts.append(m["sequence"])
                parts.append(self._random_spacer(rng))
            parts.append(self.adapters[self.order[k]])
            # Arrays are written in the forward frame only. The read simulator picks a strand per read --
            # measured at 42/57 and 16/16 on two runs -- so reverse-complementing here as well would
            # randomise orientation twice. An earlier version did, and the symptom was that arrays built
            # with p_reverse=0 still came back from skera 60 % reverse.
            yield f"a{idx:09d}", "".join(parts), members
            idx += 1

    def build_sc(self, molecules, rng, n_cdna=16, start_index=0):
        """Group 10x cDNA molecules into single-cell arrays of `n_cdna` segments.

        The nine known adapters bracket the first eight cDNAs; the remainder are concatenated directly,
        because the 10x cDNA already carries its own primers at both ends and nine adapters cannot bracket
        sixteen segments. See the note above: this reproduces skera's observed output on real single-cell
        arrays rather than claiming to know the MAS-16 adapter layout.
        """
        i, idx = 0, start_index
        n_bracketed = min(self.n_seg, n_cdna)
        while i < len(molecules):
            k = min(n_cdna, len(molecules) - i)
            if k < 1:
                break
            members = molecules[i:i + k]
            i += k
            parts = []
            for j, m in enumerate(members):
                if j < n_bracketed:
                    parts.append(self.adapters[self.order[j]])
                    parts.append(self._random_spacer(rng))
                parts.append(m["sequence"])
            # the closing adapter goes after the last bracketed cDNA's run, i.e. at the array end
            parts.append(self.adapters[self.order[n_bracketed]])
            yield f"s{idx:09d}", "".join(parts), members
            idx += 1

    @staticmethod
    def _random_spacer(rng):
        """The 1-2 bp random run skera reports as the RANDOM adapter at array ends."""
        return "".join(rng.choice("ACGT") for _ in range(rng.choice((1, 2))))


# ---------------------------------------------------------------------------
# 10x 5' v2 single-cell segment structure
#
# Recovered from the real HG002 Kinnex single-cell arrays rather than from a kit document, by orienting
# segments against the barcode whitelist and taking a column-wise consensus. 96.1 % of real segments carry
# a whitelist barcode at offset 22, and the layout that follows is exact:
#
#     [0:22]   CTACACGACGCTCTTCCGATCT    TruSeq Read 1 primer      (100 % column agreement)
#     [22:38]  16 bp cell barcode        (from the whitelist)
#     [38:48]  10 bp UMI                 (27 % agreement, i.e. random, as it should be)
#     [48:61]  TTTCTTATATGGG             10x 5' v2 TSO, found at offset 48 = 22+16+10 in 96.8 %
#     [61:..]  cDNA
#     ...      polyA                     median 29 bp, p10 26, p90 33
#     tail     GTACTCTGCGTTGATACCACTGCTT SMART primer, reverse complement
#
# The single-cell array carries SIXTEEN cDNAs, so the design specification's "16-mer" is right. The
# measurement is direct and robust: a whole array read contains a median of 16 10x TSOs (mean 15.1, p90 16)
# over 2,000 real reads, and the array averages 17,092 bp against a single-cDNA median of 952 bp.
#
# What is NOT resolved is the adapter layout. Nine adapters are identified -- the same nine as bulk, each
# appearing about once per array -- and an exhaustive 17-mer search finds no tenth: inside the blocks skera
# leaves unsegmented, every high-multiplicity motif is per-cDNA 10x structure (the R1 primer, the TSO, the
# polyA-SMART junction at ~6 copies per block), not a distinct array adapter. Nine adapters cannot bracket
# sixteen segments, so some cDNA boundaries must be adapter-free, joined directly where one cDNA's SMART
# primer meets the next one's R1 primer.
#
# `sc_array_layout` below therefore reproduces what skera demonstrably does to the real data rather than
# asserting a chemistry: given the nine adapters, skera returns about six clean single-cDNA segments per
# read plus one large multi-cDNA block, and that is what this layout produces. Replacing it with the true
# MAS-16 adapter list -- from PacBio, or from any reference segmented BAM for the single-cell kit -- is the
# one open item for this assay, and is recorded as such in docs/build-reference.md.
TENX_R1_PRIMER = "CTACACGACGCTCTTCCGATCT"
TENX_5P_TSO = "TTTCTTATATGGG"
TENX_SMART_RC = "GTACTCTGCGTTGATACCACTGCTT"
POLYA_MEDIAN, POLYA_P10, POLYA_P90 = 29, 26, 33


def tenx_segment(barcode, umi_seq, cdna, rng):
    """One full-length 10x 5' v2 cDNA molecule as Kinnex sequences it.

    Built in the orientation the real segments were oriented to, so a consumer that expects the barcode at
    offset 22 finds it there. Strand is not randomised here: the read simulator picks one per read, and
    randomising in both places was a bug the bulk builder already paid for.
    """
    # polyA is tight around its median, so a triangular draw over the measured decile range fits it better
    # than a normal and cannot go negative
    n_a = int(rng.triangular(POLYA_P10, POLYA_P90, POLYA_MEDIAN))
    return (TENX_R1_PRIMER + barcode + umi_seq + TENX_5P_TSO + cdna + "A" * n_a + TENX_SMART_RC)
