# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This repository analyzes VarChAMP dataset for the pillar project publication. The codebase is organized into three main sections, each with a standardized folder structure:

- `1_allele_collection/` - Data collection and processing
- `2_individual_assay_analyses/` - Analysis of specific assays (dual_ipa, imaging, ppi, threshold_calibration)  
- `3_integrated_assay_analyses/` - Cross-assay integration and structural analyses

Each section follows a consistent pattern:
- `1_inputs/` - Raw input data or annotations
- `2_analyses/` - Processing and analysis scripts
- `3_outputs/` - Intermediate outputs for use by other analyses

## Development Environment Setup

### Primary Environment
The main computational environment uses conda with comprehensive dependencies:

```bash
# Create and activate the main environment
conda env create -f 2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/env.yml
conda activate varchamp
```

This environment includes:
- Python 3.8 with scientific computing stack (numpy, pandas, matplotlib, seaborn, scikit-learn)
- Polars for high-performance data manipulation
- Jupyter notebook environment
- Snakemake for workflow management
- CUDA toolkit for GPU acceleration
- Specialized libraries: pycytominer, cytotable, duckdb, pyarrow

### Alternative Environment (macOS)
For macOS users, a lighter environment is available:

```bash
conda env create -f 2_individual_assay_analyses/dual_ipa/2_analyses/env-reqs.txt
```

## Key Analysis Workflows

### Imaging Analysis Pipeline
The imaging analysis uses Snakemake workflows located in `2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/`:

```bash
# Navigate to the snakemake directory
cd 2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/

# Run image preprocessing QC
python 1.image_preprocess_qc/scripts/1_calc_plate_bg.py --batch_list "batch_names" --input_dir "input_path" --output_dir "output_path" --workers 256

# Run snakemake pipeline for specific batches
snakemake --snakefile Snakefile_batch[X] --directory [working_dir] --cores 256
```

### Dual IPA Analysis
Located in `2_individual_assay_analyses/dual_ipa/2_analyses/`:
- Process flow cytometry data
- Compute allele-level statistics
- Generate perturbation scores

### PPI Analysis
Located in `2_individual_assay_analyses/ppi/2_analyses/`:
- Edgotyping score analysis
- Likelihood ratio calculations using R scripts

### Threshold Calibration
Located in `2_individual_assay_analyses/threshold_calibration/`:

```bash
# Run likelihood ratio calculations
./runLLR.sh
./runLLR_max.sh

# R-based analysis
Rscript buildLLR.R
Rscript calcLLR.R
```

## Code Architecture

### Data Processing Pipeline
1. **Raw Data Ingestion**: Scripts in `1_inputs/` directories download and prepare raw datasets
2. **Individual Assay Processing**: Each assay type has specialized processing modules
3. **Integration**: Cross-assay analysis combines results from multiple assays

### Key Python Modules

#### Utility Functions (`2_individual_assay_analyses/utils.py`)
- `plot_assay_hit_by_category()` - Visualization for hit rates by category
- `compute_aubprc()` - Adjusted precision-recall metrics
- `plot_auroc_curves()` - ROC curve analysis
- `plot_gene_level_summary()` - Gene-level hit summaries

#### Imaging Pipeline (`2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/`)
- `preprocess/` - Data cleaning and normalization modules
- `classification/` - Machine learning classification utilities
- `profiling/` - Quality control and profiling tools

## Data Standards

### File Formats
- Primary data format: TSV/CSV files
- Intermediate processing: Parquet files for performance
- Notebooks: Jupyter (.ipynb) for analysis workflows
- Configuration: JSON files for batch processing

### Naming Conventions
- Batch naming: `YYYY_MM_DD_Batch_XX` format
- Gene alleles: `gene_allele` column standard
- ClinVar categories: `clinvar_clnsig_clean` with standardized values:
  - `1_Pathogenic`, `2_Benign`, `3_Conflicting`, `4_VUS`, `5_Others`, `6_No_ClinVar`

## Running Analyses

### Jupyter Notebooks
Most analysis workflows are implemented as numbered Jupyter notebooks in `2_analyses/` directories. Run them in sequence:

```bash
# Example: imaging analysis workflow
jupyter notebook 2_individual_assay_analyses/imaging/2_analyses/2_analyze_cp_results/
# Run notebooks 0_ through 5_ in order
```

### Batch Processing
For large-scale data processing, use the shell scripts provided:

```bash
# Download raw data
./download_aws_cpg_data.sh
./download_raw_FACS.sh

# Run snakemake pipelines
./run_snakemake_pipeline.sh
```

## Git Workflow

This repository follows a fork-and-branch workflow:

```bash
# Fork the repository on GitHub, then clone your fork
git clone https://github.com/<YOUR-USERNAME>/2025_Pillar_VarChAMP.git
cd 2025_Pillar_VarChAMP

# Create feature branch
git checkout -b descriptive-branch-name

# Make changes, then commit and push
git add .
git commit -m "descriptive commit message"
git push  # Follow git's suggestions for first push

# Create pull request via GitHub interface
```

## Performance Considerations

- The imaging pipeline requires substantial computational resources (256+ cores recommended)
- GPU acceleration is available for CUDA-compatible workflows
- Large datasets use Parquet format for efficient I/O
- Snakemake provides parallel processing capabilities for batch operations