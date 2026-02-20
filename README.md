# 2025_Pillar_VarChAMP

This repository contains the code and data to reproduce **Figure 7** from the VarChAMP Pillar project publication.

**Paper:** [A Scalable Variant Effect Mapping Resource for Accelerating Comprehensive Genome Interpretation](https://www.biorxiv.org/content/10.64898/2026.02.14.705848v1.full)

## Repository Structure

```
2025_Pillar_VarChAMP/
├── 1_allele_collection/           # Allele processing from Zenodo data
├── 2_dms_bms_overlap_analyses/    # F9 imaging analysis and figures
├── 3_integrated_assay_analyses/   # Cross-assay integration figures
├── utils/                         # Shared utilities
└── env.yml                        # Conda environment
```

## Quick Start

### 1. Set up environment
```bash
conda env create -f env.yml
conda activate varchamp
```

### 2. Run analyses

**F9 Imaging Figures:**
```bash
cd 2_dms_bms_overlap_analyses/imaging/2_analyses/F9_analyses/
jupyter notebook 1_F9_visualizations.ipynb
jupyter notebook 2_F9_cell_crops.ipynb
```

**Integrated Assay Figures:**
```bash
cd 3_integrated_assay_analyses/2_analyses/
jupyter notebook 0_integrative_assay_summary.ipynb
```

## Output Figures

All manuscript figures are saved as SVG files in `3_outputs/pillar_manuscript_figures/` within each analysis directory.

## External Resources

| Resource | Description |
|----------|-------------|
| [VarChAMP Snakemake Pipeline](https://github.com/broadinstitute/2025_varchamp_snakemake/releases/tag/PillarManuscript) | Cell Painting image processing pipeline (run batches 7, 8, 11, 12, 13, 14, 15, 16 to reproduce from raw images) |
| [2025_laval_submitted](https://github.com/broadinstitute/2025_laval_submitted) | Source of integrated assay results (will be public with paper release) |

## License

See [LICENSE](LICENSE) for details.
