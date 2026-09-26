#!/usr/bin/env bash
# ============================================================
# verify_streme_ranks.sh
#
# Lists MOTIF entries in each heterodimer's repeat-masked STREME output
# (streme --dna --minw 6 --maxw 20 --thresh 0.05 on
# *_joint_peaks_filt_repeatmasked.fa), with rank and E-value. Flags entries
# whose consensus contains CACTGA or its reverse complement TCAGTG, and
# prints the first such motif per heterodimer.
#
# Usage: bash 02_motif_analysis/verify_streme_ranks.sh
# ============================================================
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
CORE_REGEX='CACTGA|TCAGTG'
SUMMARY="$(mktemp)"
trap 'rm -f "${SUMMARY}"' EXIT

echo "================================================================="
echo "STREME MOTIF RANK VERIFICATION (repeat-masked joint-peak FASTA)"
echo "================================================================="
printf "%-9s %5s %10s  %s\n" "het" "rank" "E" "consensus (first CACTGA-containing motif)" > "${SUMMARY}"
for pair in Trh_Tgo:chapter3 Sim_Tgo:chapter4 Dys_Tgo:chapter4 Sima_Tgo:chapter4; do
  het="${pair%%:*}"
  dir="${pair#*:}"
  file="${PROJECT_DIR}/data/${dir}/${het}_joint_peaks_filt_repeatmasked_streme/streme.txt"
  echo
  echo "----- ${het} -----"
  if [ ! -f "${file}" ]; then
    echo "  MISSING: ${file}"
    continue
  fi
  echo "  File: data/${dir}/${het}_joint_peaks_filt_repeatmasked_streme/streme.txt"
  echo "  Total MOTIF blocks: $(grep -c '^MOTIF ' "${file}")"
  # One line per motif: rank, consensus, E-value (from the matrix header line)
  awk -v rx="${CORE_REGEX}" -v het="${het}" -v summary="${SUMMARY}" '
    /^MOTIF / { rank++; id = $2; sub(/^[0-9]+-/, "", id); next }
    /^letter-probability matrix/ {
      e = $NF
      flag = (id ~ rx) ? "*** Motif2 candidate ***" : ""
      printf "  rank %2d  E = %-10s %-22s %s\n", rank, e, id, flag
      if (flag != "" && !done) {
        printf "%-9s %5d %10s  %s\n", het, rank, e, id >> summary
        done = 1
      }
    }' "${file}"
done
echo
echo "================================================================="
cat "${SUMMARY}"
