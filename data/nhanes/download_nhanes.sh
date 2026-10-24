#!/bin/bash
# Download the NHANES public-use files used for the thyroid analytes (cycles 2007-08, 2009-10, 2011-12).
# Run from the repository root:  bash data/nhanes/download_nhanes.sh
set -e
mkdir -p work/nht_raw
for c in "2007:E" "2009:F" "2011:G"; do
  y=${c%%:*}; s=${c##*:}
  for f in DEMO MCQ RXQ_RX THYROD; do
    curl -fsSL -o work/nht_raw/${f}_${s}.xpt "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/${y}/DataFiles/${f}_${s}.xpt"
    echo "downloaded ${f}_${s}.xpt"
  done
done
