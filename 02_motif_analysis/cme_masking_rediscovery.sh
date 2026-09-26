#!/usr/bin/env bash
# ============================================================
# cme_masking_rediscovery.sh: mask CME (ACGTG) in peak FASTA, rerun STREME to test if Motif2 is CME+AT-tract artefact.
# Input: Trh_Tgo_joint_peaks_filt_raw.fa or Trh_Tgo_joint_peaks_filt_hardmasked.fa, whichever was uploaded to MEME-ChIP.
# CME source: Long et al. 2014, PLoS One 9(1):e85518.
# Outputs: cme_positions.bed, peaks_masked.fa, streme_masked/streme.html, tomtom_vs_motif2/.
# 
set -euo pipefail

# ---- Paths (relative to PROJECT_DIR) ----
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
# Pick whichever FASTA was uploaded to MEME-ChIP. The hardmasked version
# is preferable if it was used for the original Chapter 3 motif discovery
# (because the diagnostic and the analysis must use the same sequences),
# but the raw version is acceptable if MEME-ChIP was run on it directly.
PEAKS_FA="${PROJECT_DIR}/data/chapter3/Trh_Tgo_joint_peaks_filt_raw.fa"

# Path to the Motif2 PWM in MEME format, used for the cross-heterodimer
# FIMO runs:
#   data/chapter3/STRCACTGARARAAA.meme
MOTIF2_MEME="${PROJECT_DIR}/data/chapter3/STRCACTGARARAAA.meme"

OUT_DIR="${PROJECT_DIR}/data/chapter4/motif2_diagnostic"
mkdir -p "${OUT_DIR}"

# Whether to also mask the Dys/Tgo TCGTG low-affinity variant. For a
# Trh/Tgo-specific diagnostic this should be false (Long 2014 Table 1
# lists Trh/Tgo as recognising ACGTG only).
INCLUDE_TCGTG_VARIANT="false"

# ---- Pre-flight with helpful diagnostic messages ----
if [ ! -f "${PEAKS_FA}" ]; then
  echo "MISSING: ${PEAKS_FA}"
  echo
  echo "This script expects the same peak FASTA Chapter 3 used as input"
  echo "to MEME-ChIP. Find it on your Mac:"
  echo "  find ${PROJECT_DIR} -name 'Trh_Tgo_joint_peaks*.fa' 2>/dev/null"
  echo "Then edit PEAKS_FA at the top of this script to match."
  exit 1
fi

if [ ! -f "${MOTIF2_MEME}" ]; then
  echo "MISSING: ${MOTIF2_MEME}"
  echo
  echo "This script expects the Motif2 PWM in MEME format, named"
  echo "STRCACTGARARAAA.meme (per your Chapter 4 Rmd line 2732)."
  echo "If it lives elsewhere, find it:"
  echo "  find ${PROJECT_DIR} -name 'STRCACTGARARAAA.meme' 2>/dev/null"
  echo "Then edit MOTIF2_MEME at the top of this script to match."
  exit 1
fi

for cmd in streme tomtom fimo bedtools samtools; do
  command -v "${cmd}" >/dev/null || { echo "${cmd} not in PATH"; exit 1; }
done

echo "[$(date +%H:%M:%S)] Starting Motif2 artefactual diagnostic..."
echo "[$(date +%H:%M:%S)] Input FASTA: ${PEAKS_FA}"
echo "[$(date +%H:%M:%S)] CME masked: ACGTG ($([ "${INCLUDE_TCGTG_VARIANT}" = "true" ] && echo "+ TCGTG variant" || echo "canonical only"))"
echo "[$(date +%H:%M:%S)] Output dir: ${OUT_DIR}"
echo

# ---- Step 1: Build the 5-bp CME PWM ----
cat > "${OUT_DIR}/cme.meme" << 'EOF'
MEME version 5
ALPHABET= ACGT
strands: + -
Background letter frequencies
A 0.29 C 0.21 G 0.21 T 0.29

MOTIF CME ACGTG
letter-probability matrix: alength= 4 w= 5 nsites= 100 E= 0
0.97  0.01  0.01  0.01
0.01  0.97  0.01  0.01
0.01  0.01  0.97  0.01
0.01  0.01  0.01  0.97
0.01  0.01  0.97  0.01
EOF

if [ "${INCLUDE_TCGTG_VARIANT}" = "true" ]; then
  cat >> "${OUT_DIR}/cme.meme" << 'EOF'

MOTIF CME_var TCGTG
letter-probability matrix: alength= 4 w= 5 nsites= 100 E= 0
0.01  0.01  0.01  0.97
0.01  0.97  0.01  0.01
0.01  0.01  0.97  0.01
0.01  0.01  0.01  0.97
0.01  0.01  0.97  0.01
EOF
fi

# ---- Step 2: FIMO-scan for CME instances ----
echo "[$(date +%H:%M:%S)] Scanning peaks for CME instances..."
fimo \
  --thresh 1e-3 \
  --oc "${OUT_DIR}/fimo_cme" \
  "${OUT_DIR}/cme.meme" \
  "${PEAKS_FA}"

CME_HITS="${OUT_DIR}/fimo_cme/fimo.tsv"
[ -s "${CME_HITS}" ] || { echo "FIMO returned empty output; aborting"; exit 1; }
N_HITS=$(grep -c -v "^#\|^motif_id\|^$" "${CME_HITS}" || true)
echo "[$(date +%H:%M:%S)] FIMO returned ${N_HITS} CME hits across both strands."

if [ "${N_HITS}" -lt 50 ]; then
  echo "WARNING: Fewer than 50 CME hits detected. Either the peak FASTA"
  echo "         is unusually small or CME is unexpectedly rare."
  echo "         Inspect ${CME_HITS} before trusting downstream result."
fi

# ---- Step 3: Convert FIMO hits to a BED for masking ----
awk -F'\t' 'NR>1 && $1!="" && $1 !~ /^#/ {print $3"\t"($4-1)"\t"$5}' \
    "${CME_HITS}" | sort -k1,1 -k2,2n > "${OUT_DIR}/cme_positions.bed"

N_BED_RECORDS=$(wc -l < "${OUT_DIR}/cme_positions.bed")
echo "[$(date +%H:%M:%S)] Wrote ${N_BED_RECORDS} CME positions to BED."

# ---- Step 4: Mask CME positions in the peak FASTA ----
samtools faidx "${PEAKS_FA}"

bedtools maskfasta \
  -fi "${PEAKS_FA}" \
  -bed "${OUT_DIR}/cme_positions.bed" \
  -fo "${OUT_DIR}/peaks_masked.fa" \
  -mc N

N_BEFORE=$(grep -v "^>" "${PEAKS_FA}"               | tr -cd 'N' | wc -c)
N_AFTER=$( grep -v "^>" "${OUT_DIR}/peaks_masked.fa" | tr -cd 'N' | wc -c)
N_MASKED=$((N_AFTER - N_BEFORE))
EXPECTED_BP=$((N_HITS * 5))
echo "[$(date +%H:%M:%S)] Masked ${N_MASKED} bp (expected ~${EXPECTED_BP} = ${N_HITS} hits x 5 bp; overlap-merged when adjacent)"

# ---- Step 5: Re-run STREME on the masked FASTA ----
echo "[$(date +%H:%M:%S)] Running STREME on CME-masked sequences..."
streme \
  --p     "${OUT_DIR}/peaks_masked.fa" \
  --oc    "${OUT_DIR}/streme_masked" \
  --dna \
  --nmotifs 5 \
  --minw 6 \
  --maxw 20 \
  --thresh 0.05

# ---- Step 6: TOMTOM-compare the new motifs against the original Motif2 ----
echo "[$(date +%H:%M:%S)] TOMTOM-comparing masked motifs to original Motif2..."
tomtom \
  -oc "${OUT_DIR}/tomtom_vs_motif2" \
  -thresh 0.05 \
  "${OUT_DIR}/streme_masked/streme.txt" \
  "${MOTIF2_MEME}"

# ---- Interpretation guide ----
echo
echo "=== INTERPRETATION OF RESULTS ==="
echo
echo "Inspect:"
echo "  open ${OUT_DIR}/streme_masked/streme.html"
echo "  open ${OUT_DIR}/tomtom_vs_motif2/tomtom.html"
echo
echo "CASE A — Motif2 IS recovered from CME-masked sequences:"
echo "  * STREME's top motif from the masked FASTA matches Motif2 in"
echo "    TOMTOM (q < 0.05) with similar length and AT-content."
echo "  * INTERPRETATION: Motif2 is biologically distinct from the CME."
echo "    Chapter 4's specificity claim stands without qualification."
echo
echo "CASE B — Motif2 is NOT recovered, only short AT-rich fragments:"
echo "  * STREME returns only short, low-information AT-rich motifs."
echo "  * TOMTOM does not match these to Motif2."
echo "  * INTERPRETATION: Motif2 is a CME-flanking-context artefact."
echo "    Chapter 4 needs reframing as 'enrichment of CME plus AT-rich"
echo "    flanking context' rather than a discrete novel motif."
echo
echo "CASE C — A truncated/shorter Motif2 is recovered:"
echo "  * The motif resembles Motif2 but is shorter or has lower"
echo "    information content at flanking positions."
echo "  * INTERPRETATION: Motif2 is genuinely distinct from the CME but"
echo "    has a substantial CME-overlap component. Chapter 4 needs a"
echo "    footnote acknowledging this; the central claim survives."
echo
echo "[$(date +%H:%M:%S)] Diagnostic complete."
