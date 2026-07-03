#!/usr/bin/env bash

# bcapa.sh - Bacterial Chromosome And Plasmid Assembler
#
# Single-sample hybrid Autocycler assembly restricted to two long-read
# assemblers: Flye (whole-genome) and Plassembler (plasmid-focused). Runs
# both assemblers on several subsampled read sets, then feeds all resulting
# contigs through Autocycler's compress/cluster/trim/resolve/combine steps to
# produce one consensus assembly.
#
# Assembler precedence (see the weight-tagging block below):
#   - Chromosome: Flye is the only assembler that outputs a chromosome-sized
#     contig here (Plassembler's helper output is plasmid-only), so Flye is
#     the chromosome's consensus authority by default.
#   - Plasmids: Plassembler contigs are given a higher Autocycler consensus
#     weight than Flye's, so Plassembler's plasmid calls take precedence over
#     Flye's whenever both assemblers reconstruct the same plasmid.
#
# Usage:
#   bcapa.sh <read_fastq> <threads> <jobs> [read_type]
#
#   read_fastq  Long-read FASTQ file for one sample (gzip OK).
#   threads     Threads given to each individual assembly job (capped at 128).
#   jobs        Number of assembly jobs to run simultaneously via GNU parallel
#               (8 jobs are queued: flye+plassembler x 4 subsamples, so jobs
#               beyond 8 have no additional effect).
#   read_type   One of ont_r9, ont_r10, pacbio_clr, pacbio_hifi. Default: ont_r10.
#
# Requirements:
#   - Run inside a directory containing (or reachable to write) subsampled_reads/,
#     assemblies/, and autocycler_out/ - the script writes all output relative
#     to the current working directory.
#   - conda env "autocycler" activated (autocycler, flye, plassembler, GNU
#     parallel all on PATH). See https://github.com/rrwick/Autocycler.
#
# Output:
#   autocycler_out/consensus_assembly.fasta (+ .gfa) - the final assembly.
#   autocycler_out/consensus_assembly.yaml - per-assembly stats (read by
#   `autocycler table`, as used in bcapa_pipeline.sh).
#
# Based on Ryan Wick's Autocycler example assembly script, trimmed down to
# just the flye/plassembler assembler pair and given asymmetric consensus
# weights per assembler.
#
# Copyright 2025 Ryan Wick (rrwick@gmail.com)
# Licensed under the GNU General Public License v3.
# See https://www.gnu.org/licenses/gpl-3.0.html.
# modified by Joachim Hegstad 03.11.26

# Ensure script exits on error.
set -e

#run conda init
#conda init
#use conda environment autocycler
#source activate autocycler

# Get arguments.
reads=$1                 # input reads FASTQ
threads=$2               # threads per job
jobs=$3                  # number of simultaneous jobs
read_type=${4:-ont_r10}  # read type (default = ont_r10)

# Input assembly jobs that exceed this time limit will be killed
max_time="8h"

# Validate input parameters
if [[ -z "$reads" || -z "$threads" || -z "$jobs" ]]; then
    echo "Usage: $0 <read_fastq> <threads> <jobs> [read_type]" 1>&2
    exit 1
fi
if [[ ! -f "$reads" ]]; then
    echo "Error: Input file '$reads' does not exist." 1>&2
    exit 1
fi
if (( threads > 128 )); then threads=128; fi  # Flye won't work with more than 128 threads
case $read_type in
    ont_r9|ont_r10|pacbio_clr|pacbio_hifi) ;;
    *) echo "Error: read_type must be ont_r9, ont_r10, pacbio_clr or pacbio_hifi" 1>&2; exit 1 ;;
esac

genome_size=$(autocycler helper genome_size --reads "$reads" --threads "$threads")

# Step 1: subsample the long-read set into multiple files
autocycler subsample --reads "$reads" --out_dir subsampled_reads --genome_size "$genome_size" 2>> autocycler.stderr

# Step 2: assemble each subsampled file
mkdir -p assemblies
rm -f assemblies/jobs.txt
for assembler in flye plassembler; do
    for i in 01 02 03 04; do
        echo "autocycler helper $assembler --reads subsampled_reads/sample_$i.fastq --out_prefix assemblies/${assembler}_$i --threads $threads --genome_size $genome_size --read_type $read_type" --min_depth_rel 0.1 >> assemblies/jobs.txt
    done
done
set +e
nice -n 19 parallel --jobs "$jobs" --joblog assemblies/joblog.tsv --results assemblies/logs --timeout "$max_time" < assemblies/jobs.txt
set -e

# Give all Plassembler plasmid contigs (circular or linear) extra clustering
# weight, so a plasmid found only by Plassembler still passes QC, and a
# consensus weight higher than Flye's, so Plassembler takes precedence for
# plasmids. Applied unconditionally per contig (not just circular ones): Plassembler's
# own header format is inconsistent about "circular=True" vs "circular=true",
# and linear contigs carry no circular= tag at all.
shopt -s nullglob
for f in assemblies/plassembler*.fasta; do
    sed -i 's/^>.*$/& Autocycler_cluster_weight=3 Autocycler_consensus_weight=3/' "$f"
done

# Give Flye (and Canu) contigs a consensus weight so Flye is the chromosome
# authority. Plassembler only outputs plasmid contigs, so Flye is the sole
# source of the chromosome and wins it by default; its weight of 2 here is
# below Plassembler's 3 above, so Plassembler still wins on plasmids.
for f in assemblies/canu*.fasta assemblies/flye*.fasta; do
    sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' "$f"
done
shopt -u nullglob

# Remove the subsampled reads to save space
rm subsampled_reads/*.fastq

# Step 3: compress the input assemblies into a unitig graph
autocycler compress -i assemblies -a autocycler_out 2>> autocycler.stderr

# Step 4: cluster the input contigs into putative genomic sequences
autocycler cluster -a autocycler_out 2>> autocycler.stderr

# Steps 5 and 6: trim and resolve each QC-pass cluster
for c in autocycler_out/clustering/qc_pass/cluster_*; do
    autocycler trim -c "$c" 2>> autocycler.stderr
    autocycler resolve -c "$c" 2>> autocycler.stderr
done

# Step 7: combine resolved clusters into a final assembly
autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa 2>> autocycler.stderr

#go back to previous conda environment
#conda deactivate
