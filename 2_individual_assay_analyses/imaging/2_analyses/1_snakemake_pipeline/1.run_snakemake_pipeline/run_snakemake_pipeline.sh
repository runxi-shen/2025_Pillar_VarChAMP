#!/bin/bash

## Run Image QC Pipeline
cd /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/1.image_preprocess_qc/scripts

# python 1_calc_plate_bg.py --batch_list "2025_01_27_Batch_13,2025_01_28_Batch_14" --input_dir "/home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/1.image_preprocess_qc/inputs/cpg_imgs" --output_dir "/home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/1.image_preprocess_qc/outputs/plate_bg_summary" --workers 256
# python 1_calc_plate_bg.py --batch_list "2025_06_10_Batch_18,2025_06_10_Batch_19" --input_dir "/home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/1.image_preprocess_qc/inputs/cpg_imgs" --output_dir "/home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/1.image_preprocess_qc/outputs/plate_bg_summary" --workers 256


## Run Snakemake Pipeline
cd /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/

## Run Batch 13-14
cp ../../1.run_snakemake_pipeline/snakemake_files/Snakefile_batch13 .
cp ../../1.run_snakemake_pipeline/snakemake_files/Snakefile_batch14 .

## Run batch 14
snakemake \
    --snakefile Snakefile_batch13 \
    --directory /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/ \
    --cores 256 &> /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/3_outputs/1_snakemake_pipeline/2.sm_pipeline_outputs/snakemake_logs/snakemake_batch13.log

# # Run batch 14
snakemake \
    --snakefile Snakefile_batch14 \
    --cores 256 &> /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/3_outputs/1_snakemake_pipeline/2.sm_pipeline_outputs/snakemake_logs/snakemake_batch14.log


## Run Batch 18-19
# cp snakemake_files/Snakefile_batch18 /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/
# cp snakemake_files/Snakefile_batch19 /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/

## Run batch 18
# snakemake \
#     --snakefile Snakefile_batch18 \
#     --directory /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/2_analyses/1_snakemake_pipeline/2025_varchamp_snakemake/2.snakemake_pipeline/ \
#     --cores 256 &> /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/3_outputs/1_snakemake_pipeline/2.sm_pipeline_outputs/snakemake_logs/snakemake_batch18.log

# # Run batch 18
# snakemake \
#     --snakefile Snakefile_batch19 \
#     --cores 256 &> /home/shenrunx/igvf/varchamp/2025_Pillar_VarChAMP/2_individual_assay_analyses/imaging/3_outputs/1_snakemake_pipeline/2.sm_pipeline_outputs/snakemake_logs/snakemake_batch19.log