# External Snakemake Pipeline

The Cell Painting imaging analysis in this repository uses outputs from the VarChAMP snakemake pipeline. To reproduce the full analysis from raw images, use the public release below.

## Public Release

**Repository**: [broadinstitute/2025_varchamp_snakemake](https://github.com/broadinstitute/2025_varchamp_snakemake)

**Release**: [PillarManuscript](https://github.com/broadinstitute/2025_varchamp_snakemake/releases/tag/PillarManuscript)

## Reproducing Results

1. Download the release:
   ```bash
   wget https://github.com/broadinstitute/2025_varchamp_snakemake/archive/refs/tags/PillarManuscript.tar.gz
   tar -xzf PillarManuscript.tar.gz
   cd 2025_varchamp_snakemake-PillarManuscript
   ```

2. Follow the setup instructions in the release README

3. Run the pipeline for these batches:
   - Batch 7, 8, 11, 12, 13, 14, 15, 16

## What This Analysis Uses

The F9 visualization notebooks in `../2_analyses/F9_analyses/` use pre-extracted parquet files in `../3_outputs/data/F9/`. These were generated from the pipeline outputs using `0_extract_f9_data.py`.

If you want to regenerate these files from scratch:
1. Run the snakemake pipeline (batches above)
2. Update paths in `0_extract_f9_data.py` to point to your pipeline outputs
3. Run `python 0_extract_f9_data.py`
