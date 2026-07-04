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
# Optional hairpin/topology QC (linear-plasmid-hairpin-tools): if
# HAIRPIN_TOOLS_DIR (default: ~/linear-plasmid-hairpin-tools) contains
# autocycler_dotplot_classify.py, find_hairpins.py and add_hairpin_edges.py,
# this script also runs, per sample:
#   - autocycler_dotplot_classify.py --consensus  (classify each replicon's
#     topology - circular/linear/hairpin - from the final consensus assembly)
#   - find_hairpins.py --extract-hairpins  (flag hairpin-forming contigs in
#     the raw Flye/Plassembler assemblies)
#   - add_hairpin_edges.py --overlap  (annotate the raw GFAs with hairpin
#     links so they're visible in Bandage / other topology tools)
# If HAIRPIN_TOOLS_DIR isn't found, this step is skipped with a warning.
#
# Finally, a Bandage image of each sample's consensus assembly graph is
# rendered (Bandage on PATH required; see https://github.com/rrwick/Bandage).
# Skipped with a warning if Bandage isn't found.
#
# Output:
#   ./BCAPA_OUT/<sample>_bcapa/autocycler_out/  - per-sample assembly, as
#     produced by bcapa.sh; consensus_topology.txt/.tsv are added here by the
#     hairpin/topology QC step, if it runs.
#   ./BCAPA_OUT/<sample>_bcapa/assemblies/  - raw per-assembler output;
#     hairpin_report.*, hairpin_contigs.fasta, *.hairpins.gfa,
#     hairpin_edges.tsv and hairpin_summary.tsv are added here by the
#     hairpin/topology QC step, if it runs.
#   ./BCAPA_OUT/metrics.tsv - one row per sample of Autocycler's summary
#     stats (see `autocycler table --help` for the field list).
#   ./BCAPA_OUT/<sample>_consensus_graph.png - Bandage rendering of that
#     sample's consensus_assembly.gfa.
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

# Hairpin/topology QC (see header comment above); skipped if the tools aren't found.
HAIRPIN_TOOLS_DIR="${HAIRPIN_TOOLS_DIR:-$HOME/linear-plasmid-hairpin-tools}"
if [[ -f "$HAIRPIN_TOOLS_DIR/autocycler_dotplot_classify.py" ]]; then
    for sample in *_bcapa; do
        # 1. Classify each replicon's topology in the final consensus assembly
        python3 "$HAIRPIN_TOOLS_DIR/autocycler_dotplot_classify.py" \
            --search-dir "$sample/autocycler_out" --consensus

        # 2. See which assemblers expose the hairpin on the raw outputs
        python3 "$HAIRPIN_TOOLS_DIR/find_hairpins.py" \
            --dir "$sample/assemblies" --extract-hairpins

        # 3. Annotate the raw GFAs so the hairpins show up in Bandage / topology tools
        python3 "$HAIRPIN_TOOLS_DIR/add_hairpin_edges.py" "$sample"/assemblies/*.gfa --overlap
    done
else
    echo "Warning: HAIRPIN_TOOLS_DIR ('$HAIRPIN_TOOLS_DIR') not found; skipping hairpin/topology QC." 1>&2
fi

# Render a Bandage image of each sample's consensus assembly graph.
if command -v Bandage >/dev/null 2>&1; then
    for sample in *_bcapa; do
        name="${sample%_bcapa}"
        QT_QPA_PLATFORM=offscreen Bandage image "$sample/autocycler_out/consensus_assembly.gfa" "${name}_consensus_graph.png"
    done
else
    echo "Warning: Bandage not found on PATH; skipping consensus graph images." 1>&2
fi
