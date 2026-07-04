# BCAPA — Bacterial Chromosome And Plasmid Assembler

BCAPA is a thin wrapper around [Autocycler](https://github.com/rrwick/Autocycler)
that restricts the assembler ensemble to two long-read assemblers:

- **Flye** — whole-genome assembly, used as the consensus authority for the
  bacterial **chromosome**.
- **Plassembler** — plasmid-focused assembly, used as the consensus authority
  for **plasmids**.

Both assemblers are run on several subsampled read sets, and their contigs
are combined through Autocycler's compress/cluster/trim/resolve/combine
pipeline into a single consensus assembly per sample.

## Why two assemblers, each restricted to its strength

Flye reliably assembles the chromosome but not always every plasmid;
Plassembler is built specifically for plasmid recovery and is usually more
accurate there, but doesn't assemble the chromosome at all (its Autocycler
helper output is plasmid-only). Rather than let both assemblers compete for
the same replicon in Autocycler's consensus vote, `bcapa.sh` gives each
assembler exclusive ownership of what it's good at:

- **Chromosome**: after Flye assembles each subsampled read set, `bcapa.sh`
  keeps only Flye's single largest contig and discards the rest (Flye's
  smaller contigs are typically incomplete or spurious plasmid attempts).
  Only Flye's chromosome contig reaches Autocycler.
- **Plasmids**: Plassembler's output is plasmid-only, so once Flye's plasmid
  attempts are discarded, Plassembler is the sole plasmid contributor to
  every plasmid cluster.

Because Flye and Plassembler no longer land in the same Autocycler cluster,
no consensus-weight tie-break between them is needed. Plassembler contigs
still get an
[`Autocycler_cluster_weight=3`](https://github.com/rrwick/Autocycler/wiki/Influencing-Autocycler-via-contig-headers)
tag, applied to every Plassembler contig (circular or linear), so a plasmid
found in only 1-2 of the 4 subsamples still clears Autocycler's QC threshold.

The chromosome filter is a direct FASTA filter (keep the longest contig per
`flye_*.fasta`), not `autocycler clean`/`gfa2fasta`: the `.gfa` that
Autocycler's Flye helper produces is Flye's own raw assembly graph (segments
are graph edges, not final contigs), not an Autocycler unitig graph, so it
doesn't cleanly map onto "remove this contig".

## Scripts

### `bcapa.sh`

Runs the full single-sample assembly: subsample reads, assemble with Flye and
Plassembler, restrict Flye to its chromosome contig and tag Plassembler's
contigs as described above, then compress/cluster/trim/resolve/combine via
Autocycler.

```
bcapa.sh <read_fastq> <threads> <jobs> [read_type]
```

| Argument | Description |
|---|---|
| `read_fastq` | Long-read FASTQ file for one sample (gzip OK). |
| `threads` | Threads given to each individual assembly job (capped at 128). |
| `jobs` | Number of assembly jobs run simultaneously via GNU `parallel` (8 jobs are queued: flye + plassembler x 4 subsamples). |
| `read_type` | One of `ont_r9`, `ont_r10`, `pacbio_clr`, `pacbio_hifi`. Default: `ont_r10`. |

Run it from a directory you're happy to fill with `subsampled_reads/`,
`assemblies/`, and `autocycler_out/` — all output is written relative to the
current working directory. The final assembly is
`autocycler_out/consensus_assembly.fasta` (+ `.gfa`).

### `bcapa_pipeline.sh`

Batch wrapper: runs `bcapa.sh` on every `*.fastq.gz` file in the current
directory, one sample per subdirectory, then collects per-sample stats into
one `metrics.tsv`.

```
cd <directory containing one *.fastq.gz per sample>
bcapa_pipeline.sh
```

Threads (8) and simultaneous jobs (4) are hardcoded in the `bcapa.sh "$i" 8 4`
call; edit the script to change them. Output is collected under
`BCAPA_OUT/<sample>_bcapa/` per sample, with a combined `BCAPA_OUT/metrics.tsv`
(via `autocycler table`).

If `HAIRPIN_TOOLS_DIR` (default: `~/linear-plasmid-hairpin-tools`) contains
[linear-plasmid-hairpin-tools](https://github.com/JHegstad/linear-plasmid-hairpin-tools)'s
`autocycler_dotplot_classify.py`, `find_hairpins.py` and
`add_hairpin_edges.py`, `bcapa_pipeline.sh` also runs a hairpin/topology QC
pass per sample after `metrics.tsv` is built:

- `autocycler_dotplot_classify.py --consensus` classifies each replicon in
  the final consensus assembly as circular/linear/fragmented, writing
  `consensus_topology.txt`/`.tsv` into that sample's `autocycler_out/`.
- `find_hairpins.py --extract-hairpins` flags which raw Flye/Plassembler
  contigs form a terminal hairpin (a common artifact of linear plasmids),
  writing `hairpin_report.*` and `hairpin_contigs.fasta` into that sample's
  `assemblies/`.
- `add_hairpin_edges.py --overlap` annotates the raw per-assembler GFAs with
  the detected hairpin links, writing `*.hairpins.gfa` (+ `hairpin_edges.tsv`,
  `hairpin_summary.tsv`) into that sample's `assemblies/`, so the hairpins are
  visible in Bandage.

If `HAIRPIN_TOOLS_DIR` isn't found, this step is skipped with a warning
rather than failing the pipeline.

## Requirements

- [Autocycler](https://github.com/rrwick/Autocycler) v0.6+
- Flye
- Plassembler
- GNU `parallel`

All of the above on `PATH` — developed against a conda environment named
`autocycler`.

## License

`bcapa.sh` is derived from Ryan Wick's Autocycler example assembly script and
is licensed under the GPLv3 (see [LICENSE](LICENSE)), consistent with the
original.
