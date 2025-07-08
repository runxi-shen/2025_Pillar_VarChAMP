#!/bin/bash

## Download the necessary data and files for imaging analyses pipeline from AWS Cell-Painting Gallery (CPG)

## aws cpg paths, publicly available
AWS_IMG_PATH="s3://cellpainting-gallery/cpg0020-varchamp/broad/images"
AWS_WORKSPACE_PATH="s3://cellpainting-gallery/cpg0020-varchamp/broad/workspace"
AWS_ANALYSIS_PATH="s3://cellpainting-gallery/cpg0020-varchamp/broad/workspace/analysis"

## local paths for store the aws cpg data
CPG_IMG_PATH="cpg_imgs" ## symbolic link to /data/shenrunx/igvf/varchamp/2021_09_01_VarChAMP_imgs
SNAKEMAKE_INPUT_PATH="inputs"
# BATCHES="2025_01_28_Batch_14 2025_01_27_Batch_13"
BATCHES="2025_06_10_Batch_18 2025_06_10_Batch_19"

## create directories
# mkdir -p $CPG_IMG_PATH
# mkdir -p $SNAKEMAKE_INPUT_PATH/single_cell_profiles
mkdir -p $SNAKEMAKE_INPUT_PATH/metadata/platemaps

for batch_id in $BATCHES;
do
    ## download the raw img data
    # aws s3 sync --no-sign-request "$AWS_IMG_PATH/$batch_id/images" $CPG_IMG_PATH/$batch_id/images
    aws s3 sync --no-sign-request \
        "$AWS_ANALYSIS_PATH/$batch_id" \
        "$CPG_IMG_PATH/$batch_id/analysis" \
        --exclude "*" \
        --include "**/Cells.csv" \
        --include "**/Cytoplasm.csv" \
        --include "**/Nuclei.csv" \
        --include "**/Image.csv" \
        # --recursive
        # --dry-run

    ## download the CellProfiler single-cell profiles and metadata for snakemake analysis pipeline
    # aws s3 sync --no-sign-request --exclude "*.csv" "$AWS_WORKSPACE_PATH/backend/$batch_id" $SNAKEMAKE_INPUT_PATH/single_cell_profiles/$batch_id
    # aws s3 sync --no-sign-request "$AWS_WORKSPACE_PATH/metadata/platemaps/$batch_id" $SNAKEMAKE_INPUT_PATH/metadata/platemaps/$batch_id
done