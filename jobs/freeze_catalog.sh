#!/bin/bash
# Record a sha256 for every input the designed catalog depends on and every file it produced, into
# docs/catalog-freeze.txt.
#
# The point is falsifiability: a rerun of the designer against the same design.yaml, resources and germline
# VCFs must reproduce the output checksums exactly. If it does not, something not listed here changed --
# an annotation build, a netMHCpan version, a resource edited in place -- and the catalog is not the one the
# release was built from.
#
# Run it after any designer rerun that is meant to be kept.
set -euo pipefail
source "${IGI_ENV:?export IGI_ENV=/path/to/jobs/env.sh first}"
cd "$IGI_REPO"
OUT=docs/catalog-freeze.txt
P=$IGI_PATHS

vcfs=$(python3 - "$P" <<'PY'
import sys, yaml
p = yaml.safe_load(open(sys.argv[1]))
seen, out = set(), []
for key in ("germline_vcf", "germline_vcf_baseline"):
    for _k, v in sorted((p.get(key) or {}).items()):
        if v not in seen:            # the dataset and its baseline alias point at the same file
            seen.add(v); out.append(v)
print("\n".join(out))
PY
)

{
  echo "# IGI-SYN-SEQ catalog freeze"
  echo "#"
  echo "# Every input the designed catalog depends on, and every file it produced, with a sha256. A rerun of"
  echo "# the designer against the same design.yaml, the same resources and the same germline VCFs must"
  echo "# reproduce the output checksums exactly; if it does not, something not listed here has changed."
  echo "#"
  echo "# Regenerate with: jobs/freeze_catalog.sh"
  echo "# Frozen: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "# Catalog commit: $(git rev-parse HEAD)"
  # The image identity and the binary, not the whole command: the command carries a site-specific path and
  # no tracked file may contain one. What the freeze needs to pin is the netMHCpan VERSION.
  echo "# netMHCpan image:  $(grep -oE 'netmhcpan_cmd.*' "$P" | grep -oE '[^/]+\.(img|sif)' | head -1)"
  echo "# netMHCpan binary: $(grep -oE 'netmhcpan_cmd.*' "$P" | grep -oE '/netMHCpan[^ ]*/netMHCpan' | head -1)"
  echo
  echo "## Germline VCFs (inputs, not tracked)"
  while read -r f; do
    [ -n "$f" ] || continue
    [ -s "$f" ] || { echo "MISSING  $f"; continue; }
    printf "%s  %s\n" "$(sha256sum "$f" | cut -d' ' -f1)" "$(basename "$f")"
  done <<< "$vcfs"
  echo
  echo "## Catalog design parameters, resources and outputs (tracked)"
  ( cd catalog && sha256sum output/*.tsv output/*.json resources/* design.yaml )
} > "$OUT"
echo "wrote $OUT ($(grep -cvE '^(#|$|##)' "$OUT") checksums)"
