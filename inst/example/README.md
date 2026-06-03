# cellgpsr Example Templates

These scripts are repository-level templates for manuscript-oriented analyses.
They are not included in the installed R package build.

The examples intentionally do not bundle large public datasets. Set the
environment variables below, or edit the script-local path variables, before
running a template:

- `CELLGPSR_FIBROSIS_SEURAT_RDS`: fibrosis Seurat RDS used by the fibrosis
  examples.
- `CELLGPSR_FIBROSIS_REGION_CSV`: fibrosis region-annotation CSV.
- `CELLGPSR_NPC_TRANSCRIPTS_CSV`: NPC transcript CSV.
- `CELLGPSR_NPC_CELLS_CSV`: NPC cell-type CSV.
- `CELLGPSR_LN_CELLS_CSV`: 10x lymph-node cells CSV or CSV.GZ.
- `CELLGPSR_LN_CELLTYPES_CSV`: 10x lymph-node cell-type annotation CSV.
- `CELLGPSR_LN_TRANSCRIPTS_PARQUET`: 10x lymph-node transcripts parquet.
- `CELLGPSR_TBC_RESULT_DIR`: transcript-by-cell result directory.
- `CELLGPSR_MORPHOLOGY_IMAGE`: morphology image used by the TRU overlay
  template.
- `CELLGPSR_OUTPUT_DIR`: output directory for generated plots and CSV files.
