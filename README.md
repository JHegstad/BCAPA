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

## Why two assemblers with asymmetric weights

Autocycler builds its consensus sequence per replicon (chromosome, each
plasmid) from a weighted vote across all contigs that cluster together. Flye
is a reliable whole-genome assembler but can struggle to fully resolve small
or low-copy plasmids; Plassembler is built specifically for plasmid recovery
and is usually more accurate there, but doesn't assemble the chromosome at
all (its Autocycler helper output is plasmid-only).

`bcapa.sh` tags each assembler's contigs with an
[`Autocycler_consensus_weight`](https://github.com/rrwick/Autocycler/wiki/Influencing-Autocycler-via-contig-headers)
so that:

- Flye contigs get `Autocycler_consensus_weight=2`. Since Plassembler never
  outputs a chromosome-sized contig, Flye is the chromosome's only
  contributor and wins it by default.
- Plassembler contigs get `Autocycler_consensus_weight=3` (plus
  `Autocycler_cluster_weight=3`, so a plasmid found only by Plassembler still
  passes Autocycler's QC threshold). This weight beats Flye's, so
  Plassembler's plasmid sequence takes precedence whenever both assemblers
  reconstruct the same plasmid.

The weight is applied to every Plassembler contig regardless of whether it's
circular or linear — Plassembler's own output is inconsistent about how it
tags circular contigs (`circular=True` vs `circular=true`, and linear contigs
carry no `circular=` tag at all), so matching on that text isn't reliable.

## Scripts

### `bcapa.sh`

Runs the full single-sample assembly: subsample reads, assemble with Flye and
Plassembler, apply the weight tags above, then compress/cluster/trim/resolve/
combine via Autocycler.

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
