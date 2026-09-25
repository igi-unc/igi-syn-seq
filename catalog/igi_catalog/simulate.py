"""Stage 4: read simulation.

Reads are produced by mixing clone genomes and transcriptomes in the proportions the tumour model
implies, then handing the mixture to a simulator:

  * WES and WGS: fragments are drawn from clone genome sequences, weighted by the exclusive fraction of
    each clone and by tumour purity, then sequenced with ART using a profile trained on real data;
  * bulk RNA: molecules are drawn per transcript from the per-clone TPM table, fragmented, and sequenced
    the same way.

The mixing weights are computed here rather than inside the simulator, so the expected coverage of any
allele can be derived from the truth table without running anything.
"""
import os
import subprocess


def clone_weights(clones, purity):
    """Fraction of sequenced molecules coming from each clone, plus the normal fraction.

    A cell contributes in proportion to the exclusive fraction of its most-derived clone; normal cells
    contribute 1 - purity.
    """
    w = {c: purity * f for c, f in clones.excl.items() if f > 0}
    w["NORMAL"] = 1.0 - purity
    total = sum(w.values())
    return {k: v / total for k, v in w.items()}


def capture_weights(intervals, chrom, length, on_target=0.9):
    """Probability mass on captured versus off-target positions for an exome library."""
    return {"on_target": on_target, "off_target": 1.0 - on_target}


class ArtRunner:
    """Illumina read simulation through ART."""

    def __init__(self, cmd_template, workdir):
        self.cmd = cmd_template
        self.workdir = workdir
        os.makedirs(workdir, exist_ok=True)

    def paired(self, in_fasta, out_prefix, fold_coverage, read_len=150, frag_mean=350, frag_sd=60, seed=1,
               profile=None, extra=""):
        """Run ART in paired-end mode over `in_fasta`.

        Returns (r1_path, r2_path).
        """
        args = (f"-ss HS25 -i {in_fasta} -p -l {read_len} -f {fold_coverage} "
                f"-m {frag_mean} -s {frag_sd} -rs {seed} -na -o {out_prefix}")
        if profile:
            args += f" -1 {profile}1.txt -2 {profile}2.txt"
        cmd = self.cmd.format(args=args) + (" " + extra if extra else "")
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        return f"{out_prefix}1.fq", f"{out_prefix}2.fq"


def write_region_fasta(genome_seq, chrom, path, name=None, line=60):
    """Write one sequence to FASTA."""
    with open(path, "w") as fh:
        fh.write(f">{name or chrom}\n")
        for i in range(0, len(genome_seq), line):
            fh.write(genome_seq[i:i + line] + "\n")
    return path


def transcript_fasta_with_abundance(records, clone, path, min_tpm=0.01):
    """FASTA of transcripts expressed in `clone`, with TPM in the header.

    Returns [(id, tpm)] so the caller can turn TPM into molecule counts.
    """
    out = []
    with open(path, "w") as fh:
        for r in records:
            if r.get("clone") not in (clone, None):
                continue
            tpm = float(r.get("tpm") or 0)
            if tpm < min_tpm or not r.get("sequence"):
                continue
            fh.write(f">{r['id']} tpm={tpm}\n")
            seq = r["sequence"]
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
            out.append((r["id"], tpm))
    return out


def molecules_from_tpm(abundances, total_molecules):
    """Turn TPM into integer molecule counts summing to `total_molecules`."""
    total_tpm = sum(t for _i, t in abundances) or 1.0
    counts = []
    running = 0
    for i, (tid, tpm) in enumerate(abundances):
        n = int(round(total_molecules * tpm / total_tpm))
        if i == len(abundances) - 1:
            n = max(0, total_molecules - running)
        counts.append((tid, n))
        running += n
    return counts
