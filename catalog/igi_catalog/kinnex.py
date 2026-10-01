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

    @staticmethod
    def _random_spacer(rng):
        """The 1-2 bp random run skera reports as the RANDOM adapter at array ends."""
        return "".join(rng.choice("ACGT") for _ in range(rng.choice((1, 2))))
