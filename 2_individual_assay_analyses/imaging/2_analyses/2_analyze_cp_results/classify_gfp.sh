#!/bin/bash

cd ../1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/classification

python classify_gfp_filtered.py \
    --allele_file "/home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/F9_variants.csv" \
    --batch "2025_01_27_Batch_13" \
    --input_dir "../outputs/batch_profiles" \
    --output_dir "/home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/1_allele_collection/1_inputs/raw_inputs"

python classify_gfp_filtered.py \
    --allele_file "/home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/F9_variants.csv" \
    --batch "2025_01_28_Batch_14" \
    --input_dir "../outputs/batch_profiles" \
    --output_dir "/home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/1_allele_collection/1_inputs/raw_inputs"