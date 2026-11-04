# Yeast_test — power benchmark on the yeast 48x48 data

The yeast *Δsnf2* vs wild-type experiment of Schurch et al. (2016) (bartongroup/profDGE48):
48 biological replicates per condition, a clean reference for small-n behaviour next to TCGA
(TCGA_test/) and simulated data (Simulated_data/).

## Data (on the server)

`Distribution_GOF/slurm/gof_yeast48_submit.sh setup` built the datasets from the htseq counts:

| file | content |
|---|---|
| `/media/BioNAS2/TCGA_TOX_TEST/data/yeast48/yeast48_clean.rds` | 7126 genes x 86 samples (42 WT, 44 Snf2), bad replicates removed (default here) |
| `/media/BioNAS2/TCGA_TOX_TEST/data/yeast48/yeast48_all.rds` | all 96 samples, bad replicates flagged in `samples$bad_replicate` |

`scripts/make_yeast_cohorts.R` turns a dataset into one base cohort per condition in the format of
`TCGA_test/scripts/power_test.R` (`load_real_cohort()`): `cohorts_<dataset>/yeast48_wt.rds` and
`yeast48_snf2.rds` next to the data. The counts carry no gene lengths: pass a yeast GTF
(`--GTF`, union exon length per gene) or a table (`--LENGTHS`); without either every length is 1,
so TPM = CPM (recorded as `length_source`). The slurm scripts build the cohorts on first use.

## Runs

The engine is shared with TCGA (`TCGA_test/scripts/power_test.R`, `power_report.R`,
`knn_config_test.R`); nothing is copied here. Parts R (random halves of one condition) and RS
(2n random samples of one condition), with known log2 FCs injected by binomial thinning, work on
the yeast cohorts unchanged. Base cohorts: WT and Snf2 (2, where TCGA has 6; the report's
"5 of 6 cohorts" gates become "1 of 2").

```
mkdir -p logs
sbatch Yeast_test/slurm/yeast_power.sbatch     # Parts RS + R + report  -> $ROOT/yeast_power_out
sbatch Yeast_test/slurm/yeast_knn.sbatch       # kNN configuration grid -> $ROOT/yeast_knn_out
# options: --export=ALL,YEAST_DATASET=<...>/yeast48_all.rds,YEAST_GTF=<gtf>,REBUILD=1
```

Rebuild the TOX library (Apptainer `./build.sh`) before the first run if the Fortran changed.
Copy results to `Yeast_test/results/` (not versioned) for local analysis.
