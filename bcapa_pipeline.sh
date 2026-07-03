#!/bin/bash
#
# bcapa_pipeline.sh - batch wrapper for bcapa.sh
#
# Runs bcapa.sh (Flye + Plassembler Autocycler assembly) on every
# *.fastq.gz file in the current directory, one sample at a time, then
# collects all per-sample stats into a single metrics.tsv.
#
# Usage:
#   cd <directory containing one *.fastq.gz file per sample>
#   bcapa_pipeline.sh
#
#   Threads (8) and simultaneous jobs (4) are hardcoded below; edit the
#   `bcapa.sh "$i" 8 4` call to change them. Read type defaults to ont_r10
#   (bcapa.sh's default) since it is not passed through here.
#
# Requirements:
#   - bcapa.sh and autocycler on PATH (conda env "autocycler" activated).
#   - Run from a directory that only contains the *.fastq.gz inputs for this
#     batch; the script creates a <sample>_bcapa/ working directory per file.
#
# Output:
#   ./BCAPA_OUT/<sample>_bcapa/autocycler_out/  - per-sample assembly, as
#     produced by bcapa.sh.
#   ./BCAPA_OUT/metrics.tsv - one row per sample of Autocycler's summary
#     stats (see `autocycler table --help` for the field list).
#
# Joachim Hegstad 03.07.26
for i in *.fastq.gz
do
  OUT="$(basename $i .fastq.gz)_bcapa"
  mkdir $OUT
  cp $i $OUT
  cd $OUT
  bcapa.sh "$i" 8 4
  cd ..
done
# move into one folder and create metrics.tsv
mkdir BCAPA_OUT
mv *_bcapa BCAPA_OUT
cd BCAPA_OUT
autocycler table > metrics.tsv  # create the TSV header
for sample in *_bcapa; do
    autocycler table -a "$sample" -n "$sample" >> metrics.tsv  # append a TSV row
done
