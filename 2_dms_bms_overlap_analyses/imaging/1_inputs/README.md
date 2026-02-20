# Imaging Input Data

Input files for imaging analysis are located in the shared inputs directory:

```
../../3_integrated_assay_analyses/1_inputs/
```

## Required Files

| File | Location | Description |
|------|----------|-------------|
| `VarChAMP_data_supp_mat_PP.tsv` | `../../3_integrated_assay_analyses/1_inputs/VarChAMP_data_supp_mat_PP.tsv` | Publication-ready integrated assay results |

## Usage in Scripts

Scripts in `2_analyses/F9_analyses/` reference these files using relative paths:
- `0_extract_f9_data.py`: Uses `../../../../3_integrated_assay_analyses/1_inputs/VarChAMP_data_supp_mat_PP.tsv`
- `1_F9_visualizations.ipynb`: Uses `../../../../3_integrated_assay_analyses/1_inputs/VarChAMP_data_supp_mat_PP.tsv`

Input annotation files are shared across multiple assay analyses (imaging, dual IPA, PPI) to:
1. Maintain a single source of truth for variant annotations
2. Ensure consistency across analyses
3. Reduce data duplication
