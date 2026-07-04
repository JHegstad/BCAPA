# Validation: Flye chromosome-only restriction

Full `bcapa_pipeline.sh` run on 4 real ONT samples, comparing the pipeline
before and after restricting Flye to its single largest (chromosome) contig
and dropping the cross-assembler consensus-weight tie-break (see
`ce46c9a`).

| Sample | QC pass/fail (old → new) | Unitigs (old → new) | Fully resolved (old → new) | Consensus bases (old → new) |
|---|---|---|---|---|
| 51525510 | 11 / 27 → 10 / 0 | 22 → 10 | false → **true** | 3,251,340 → 3,243,097 |
| 51535247 | 7 / 5 → 7 / 0 | 42 → 7 | false → **true** | 3,224,891 → 3,216,562 |
| 51558664 | 11 / 13 → 10 / 0 | 18 → 10 | false → **true** | 3,317,202 → 3,297,481 |
| 51560653 | 8 / 12 → 8 / 1 | 11 → 8 | false → **true** | 3,230,955 → 3,209,448 |

All 4 samples went from a fragmented, unresolved consensus assembly to a
fully resolved one (one sequence per replicon: a circular chromosome plus
each plasmid), with QC-failed clusters eliminated or nearly eliminated.
Total assembly length is essentially unchanged (within ~0.7%), confirming
the cleanup removed fragmentation/duplication rather than real sequence.

Raw data: `metrics.tsv` from both runs, plus this run's per-sample
`consensus_topology.tsv`/`hairpin_report.tsv` and `<sample>_consensus_graph.png`
Bandage renderings, are kept alongside the test data
(`~/NGS/dev/testdata/BCAPA_OUT` for this run, `BCAPA_OUT_20260703` for the
prior one) rather than in this repo.
