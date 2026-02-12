#!/usr/bin/env python
"""
Extract F9 single-cell GFP features and feature importance data.

This script extracts:
1. Feature importance per variant (from classification results)
2. Single-cell GFP features (minimal set: essential biological features + important features)
3. Nuc2Cyto_Ratio_Mean computed from raw profiles
4. Z-score normalized variant-level features (variant vs reference)

Output files are saved with zstd compression for GitHub storage (~55MB total).

Usage:
    python extract_f9_data.py                    # Run full extraction
    python extract_f9_data.py --skip-features    # Only extract feature importance
    python extract_f9_data.py --skip-importance  # Only extract single-cell features
"""

import polars as pl
from pathlib import Path
from tqdm import tqdm
import argparse
import os


# Configuration
GENE = "F9"
FEAT = "GFP"
IMPORTANCE_THRESHOLD = 1e-2

# Essential biological features for F9 analysis (protein localization, expression, texture)
ESSENTIAL_BIOLOGICAL_FEATURES = [
    'Cells_Correlation_RWC_AGP_GFP',           # Correlation between AGP and GFP channels
    'Cells_Intensity_MeanIntensity_GFP',       # Mean GFP intensity in cells
    'Cytoplasm_Texture_Entropy_GFP_10_00_256', # GFP texture entropy in cytoplasm
    'Nuclei_Intensity_MeanIntensity_GFP',      # Mean GFP intensity in nuclei
    'Nuclei_Correlation_K_GFP_Mito',           # GFP-Mito correlation in nuclei
    'Nuclei_Correlation_K_GFP_AGP',            # GFP-AGP correlation in nuclei
    'Nuclei_Texture_DifferenceVariance_GFP_5_00_256',  # GFP texture variance in nuclei
    'Cytoplasm_Granularity_1_GFP',             # GFP granularity in cytoplasm
]

# Essential metadata columns for visualization
ESSENTIAL_METADATA = ['Metadata_CellID', 'Metadata_gene_allele', 'Metadata_Plate', 'Metadata_Well']

# Directories
BATCH_PROFILE_DIR = Path("/home/shenrunx/igvf/varchamp/2025_varchamp_snakemake/2.snakemake_pipeline/outputs/batch_profiles")
CLASSIFICATION_DIR = Path("/home/shenrunx/igvf/varchamp/2025_varchamp_snakemake/2.snakemake_pipeline/outputs/classification_results")
OUTPUT_DIR = Path("../../3_outputs/2_results_summary/F9")
INPUT_CLINVAR = Path("../../1_inputs/0_imaging_phenotype_summary_clinvar.tsv")

BIO_REP_BATCHES_DICT = {
    "2025_01_Batch_13-14": ("2025_01_27_Batch_13", "2025_01_28_Batch_14")
}


def load_clinvar_summary() -> pl.DataFrame:
    """Load the imaging summary with ClinVar annotations."""
    return pl.read_csv(INPUT_CLINVAR, separator="\t", infer_schema_length=100000)


def compute_nuc2cyto_ratio(batch_ids: tuple) -> pl.DataFrame:
    """Compute nucleus-to-cytoplasm GFP intensity ratio from raw profiles."""
    dfs = []
    for batch_id in batch_ids:
        profile_path = BATCH_PROFILE_DIR / batch_id / "profiles.parquet"
        df = (
            pl.scan_parquet(profile_path)
            .filter(pl.col("Metadata_symbol").str.contains(GENE))
            .with_columns(
                (pl.col("Nuclei_Intensity_IntegratedIntensity_GFP") /
                 pl.col("Cytoplasm_Intensity_IntegratedIntensity_GFP")).alias("Nuc2Cyto_Ratio_Integrated"),
                (pl.col("Nuclei_Intensity_MeanIntensity_GFP") /
                 pl.col("Cytoplasm_Intensity_MeanIntensity_GFP")).alias("Nuc2Cyto_Ratio_Mean")
            )
            .select([
                "Metadata_Plate", "Metadata_Well", "Metadata_ImageNumber", "Metadata_ObjectNumber",
                "Nuc2Cyto_Ratio_Integrated", "Nuc2Cyto_Ratio_Mean"
            ])
            .collect()
        )
        dfs.append(df)

    result = pl.concat(dfs, how="diagonal")
    result = result.with_columns(
        pl.concat_str(
            ["Metadata_Plate", "Metadata_Well", "Metadata_ImageNumber", "Metadata_ObjectNumber"],
            separator="_"
        ).alias("Metadata_CellID")
    )
    return result


def extract_feature_importance(imaging_summary: pl.DataFrame) -> pl.DataFrame:
    """Extract feature importance per variant from classification results."""
    feat_import_df = pl.DataFrame()

    print("Loading feature importance data from batches...")
    for bio_rep, bio_rep_batches in BIO_REP_BATCHES_DICT.items():
        allele_batch_df = imaging_summary.filter(pl.col("Metadata_Bio_Batch") == bio_rep)
        batch_allele_list = list(allele_batch_df["gene_variant"]) + list(allele_batch_df["symbol"].unique())

        for batch_id in bio_rep_batches:
            feat_importance_path = (
                CLASSIFICATION_DIR / batch_id /
                "profiles_tcdropped_filtered_var_mad_outlier_featselect_filtcells" /
                "feat_importance_gfp_adj.csv"
            )

            if not feat_importance_path.exists():
                print(f"  Warning: {feat_importance_path} not found, skipping")
                continue

            feat_df_b = pl.scan_csv(feat_importance_path)
            schema = feat_df_b.collect_schema()
            meta_cols = [c for c in schema.names() if c.startswith("Metadata_") or c.startswith("Group") or c == "Batch"]
            feat_cols = [c for c in schema.names() if c not in meta_cols and FEAT in c]

            feat_df_batch = (
                feat_df_b
                .with_columns([pl.col(c).cast(pl.Float64, strict=False).alias(c) for c in feat_cols])
                .with_columns(
                    pl.lit(batch_id).alias("Metadata_Batch"),
                    pl.col("Group2").str.split('_').list.slice(-2).list.join("_").alias("Metadata_gene_allele")
                )
                .filter(
                    (pl.col("Metadata_Feature_Type") == FEAT) & (~pl.col("Metadata_Control")) &
                    (pl.col("Group1").is_in(batch_allele_list) | pl.col("Metadata_gene_allele").is_in(batch_allele_list))
                )
                .select(meta_cols + feat_cols + ["Metadata_Batch", "Metadata_gene_allele"])
                .collect()
            )
            feat_import_df = pl.concat([feat_import_df, feat_df_batch], how="diagonal")

    # Aggregate feature importance per variant
    print("Aggregating feature importance per variant...")
    feat_df = pl.DataFrame()
    gene_variants = imaging_summary.filter(pl.col("Gene") == GENE)["gene_variant"].unique().to_list()

    for allele in tqdm(gene_variants, desc="Processing variants"):
        batch_list = list(feat_import_df.filter(
            pl.col("Metadata_gene_allele") == allele
        ).unique("Metadata_Batch")["Metadata_Batch"])

        allele_df = pl.DataFrame()
        for batch in batch_list:
            batch_df = feat_import_df.filter(
                (pl.col("Metadata_Batch") == batch) & (pl.col("Metadata_gene_allele") == allele)
            )
            feat_cols = [col for col in batch_df.columns if "Metadata" not in col and "Group" not in col]
            batch_df = batch_df.filter(~pl.all_horizontal(pl.col(feat_cols).is_null()))
            non_null_gfp = [col for col in feat_cols if not batch_df[col].is_null().any()]

            if not non_null_gfp:
                continue

            batch_df_col = (
                batch_df.select(non_null_gfp).mean()
                .transpose(include_header=True)
                .filter(pl.col("column_0") > IMPORTANCE_THRESHOLD)
                .sort("column_0", descending=True)
            )

            if not allele_df.is_empty():
                allele_df = allele_df.join(
                    batch_df_col.with_columns(
                        pl.col("column_0").alias(batch),
                        pl.col("column").alias("index"),
                    ).select("index", batch),
                    on="index", how="full", coalesce=True
                )
            else:
                allele_df = batch_df_col.with_columns(
                    pl.col("column_0").alias(batch),
                    pl.col("column").alias("index"),
                ).select("index", batch)

        if allele_df.is_empty():
            continue

        allele_df = allele_df.with_columns(
            pl.mean_horizontal(pl.col([col for col in allele_df.columns if "Batch" in col])).alias("feat_importance"),
            pl.col("index").alias("cp_feature")
        ).select(["cp_feature", "feat_importance"]).rename({"feat_importance": allele})

        if not feat_df.is_empty():
            feat_df = feat_df.join(allele_df, on="cp_feature", how="full", coalesce=True)
        else:
            feat_df = allele_df

    return feat_df


def get_important_features(feat_importance_df: pl.DataFrame) -> list:
    """Get features with importance > threshold for any variant."""
    variant_cols = [c for c in feat_importance_df.columns if c != 'cp_feature']
    feat_with_max = feat_importance_df.with_columns(
        pl.max_horizontal(variant_cols).alias('max_importance')
    )
    return feat_with_max.filter(pl.col('max_importance') > IMPORTANCE_THRESHOLD)['cp_feature'].to_list()


def extract_minimal_features(
    input_parquet: Path,
    nuc2cyto_df: pl.DataFrame,
    important_features: list
) -> pl.DataFrame:
    """Extract minimal feature set from single-cell profiles."""
    print(f"Reading {input_parquet}...")
    df = pl.read_parquet(input_parquet)

    # Combine essential biological features + important features + Nuc2Cyto_Ratio_Mean
    all_needed_features = list(set(
        ESSENTIAL_BIOLOGICAL_FEATURES + important_features + ['Nuc2Cyto_Ratio_Mean']
    ))

    # Filter to features that exist in the data
    feat_cols_present = [c for c in all_needed_features if c in df.columns]
    # Nuc2Cyto_Ratio_Mean comes from nuc2cyto_df, not the main df
    feat_cols_present = [c for c in feat_cols_present if c != 'Nuc2Cyto_Ratio_Mean']

    print(f"  Essential biological features: {len(ESSENTIAL_BIOLOGICAL_FEATURES)}")
    print(f"  Important features (threshold {IMPORTANCE_THRESHOLD}): {len(important_features)}")
    print(f"  Total unique features to extract: {len(feat_cols_present) + 1}")  # +1 for Nuc2Cyto_Ratio_Mean

    # Select metadata + feature columns
    metadata_cols_present = [c for c in ESSENTIAL_METADATA if c in df.columns]
    df_minimal = df.select(metadata_cols_present + feat_cols_present)

    # Join with Nuc2Cyto ratio
    df_minimal = df_minimal.join(
        nuc2cyto_df.select(["Metadata_CellID", "Nuc2Cyto_Ratio_Mean"]),
        on="Metadata_CellID",
        how="left"
    )

    # Convert to Float32 for compression
    float_cols = [c for c in df_minimal.columns if c not in metadata_cols_present]
    df_minimal = df_minimal.with_columns([
        pl.col(c).cast(pl.Float32) for c in float_cols
    ])

    return df_minimal


def compute_zscore_variant_features(
    sc_features_df: pl.DataFrame,
    imaging_summary: pl.DataFrame
) -> pl.DataFrame:
    """
    Compute z-score normalized features per variant (variant vs reference).

    For each variant, computes z-score of features relative to the reference (wild-type).
    Aggregates by well first (median), then by variant (median).

    Returns:
        DataFrame with one row per variant, z-score normalized features as columns.
    """
    print("Computing z-score normalized variant features...")

    # Get feature columns (non-metadata, GFP-related)
    feature_cols = [c for c in sc_features_df.columns
                    if "Metadata" not in c and ("GFP" in c or "Nuc2Cyto" in c)]

    # Get reference (wild-type) statistics
    ref_df = sc_features_df.filter(pl.col("Metadata_gene_allele") == GENE)
    ref_mean = ref_df.select(feature_cols).mean().to_numpy()[0].tolist()
    ref_std = ref_df.select(feature_cols).std().to_numpy()[0].tolist()

    # Get valid variants from imaging summary
    valid_variants = imaging_summary.filter(
        pl.col("Gene") == GENE
    )["gene_variant"].unique().to_list()

    # Filter to variant cells only (not reference)
    variant_df = sc_features_df.filter(
        (pl.col("Metadata_gene_allele").str.starts_with(f"{GENE}_")) &
        (pl.col("Metadata_gene_allele") != GENE) &
        (pl.col("Metadata_gene_allele").is_in(valid_variants))
    )

    # Z-score normalize features
    z_norm_sc_df = variant_df.with_columns([
        ((pl.col(col) - pl.lit(ref_mean[i])) / pl.lit(ref_std[i])).alias(col)
        for i, col in enumerate(feature_cols)
    ])

    # Aggregate by well first (median), then by variant (median)
    z_norm_variant_df = z_norm_sc_df.group_by(
        "Metadata_gene_allele", "Metadata_Plate", "Metadata_Well"
    ).agg([
        pl.col(col).median().alias(col) for col in feature_cols
    ]).group_by("Metadata_gene_allele").agg([
        pl.col(col).median().alias(col) for col in feature_cols
    ])

    # Join with imaging summary metadata
    metadata_cols = ["gene_variant", "img_pheno_score", "img_misloc_score",
                     "clinvar_clnsig_clean", "clinvar_clnsig_clean_pp_strict"]
    metadata_cols_present = [c for c in metadata_cols if c in imaging_summary.columns]

    z_norm_with_meta = z_norm_variant_df.join(
        imaging_summary.filter(pl.col("Gene") == GENE).select(metadata_cols_present).unique("gene_variant"),
        left_on="Metadata_gene_allele",
        right_on="gene_variant",
        how="left"
    )

    # Rename metadata columns with Metadata_ prefix
    rename_dict = {c: f"Metadata_{c}" for c in metadata_cols_present if c != "gene_variant"}
    z_norm_with_meta = z_norm_with_meta.rename(rename_dict)

    # Convert to Float32
    float_cols = [c for c in z_norm_with_meta.columns if "Metadata" not in c]
    z_norm_with_meta = z_norm_with_meta.with_columns([
        pl.col(c).cast(pl.Float32) for c in float_cols
    ])

    print(f"  Variants processed: {len(z_norm_with_meta)}")
    print(f"  Features: {len(feature_cols)}")

    return z_norm_with_meta


def main():
    parser = argparse.ArgumentParser(description="Extract F9 GFP features and feature importance")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Output directory")
    parser.add_argument("--skip-features", action="store_true", help="Skip feature extraction")
    parser.add_argument("--skip-importance", action="store_true", help="Skip feature importance extraction")
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load ClinVar summary
    print("Loading ClinVar summary...")
    imaging_summary = load_clinvar_summary()

    # Step 1: Extract feature importance (needed to determine which features to keep)
    importance_path = output_dir / "f9_feature_importance.parquet"

    if not args.skip_importance:
        print("\n" + "="*60)
        print("Step 1: Extracting feature importance per variant...")
        print("="*60)
        feat_importance_df = extract_feature_importance(imaging_summary)

        print(f"\nSaving feature importance to {importance_path}...")
        feat_importance_df.write_parquet(importance_path, compression="zstd", compression_level=22)
        size_mb = os.path.getsize(importance_path) / (1024 * 1024)
        print(f"  Output size: {size_mb:.2f} MB")
    else:
        print("\nLoading existing feature importance...")
        feat_importance_df = pl.read_parquet(importance_path)

    # Step 2: Extract minimal single-cell features
    if not args.skip_features:
        print("\n" + "="*60)
        print("Step 2: Extracting minimal single-cell features...")
        print("="*60)

        # Get important features from the importance file
        important_features = get_important_features(feat_importance_df)

        # Compute Nuc2Cyto ratio from raw profiles
        print("\nComputing Nuc2Cyto ratio from raw profiles...")
        batch_ids = BIO_REP_BATCHES_DICT["2025_01_Batch_13-14"]
        nuc2cyto_df = compute_nuc2cyto_ratio(batch_ids)

        # Extract minimal features
        input_parquet = output_dir / "ref_var_sc_profiles_mad_orig.parquet"
        df_minimal = extract_minimal_features(input_parquet, nuc2cyto_df, important_features)

        # Save compressed parquet
        output_path = output_dir / "f9_sc_features_minimal.parquet"
        print(f"\nSaving minimal features to {output_path}...")
        df_minimal.write_parquet(output_path, compression="zstd", compression_level=22)

        size_mb = os.path.getsize(output_path) / (1024 * 1024)
        print(f"  Shape: {df_minimal.shape}")
        print(f"  Output size: {size_mb:.2f} MB")

    # Step 3: Compute z-score normalized variant-level features
    if not args.skip_features:
        print("\n" + "="*60)
        print("Step 3: Computing z-score normalized variant features...")
        print("="*60)

        z_norm_df = compute_zscore_variant_features(df_minimal, imaging_summary)

        # Save compressed parquet
        zscore_path = output_dir / "f9_variant_zscore_features.parquet"
        print(f"\nSaving z-score features to {zscore_path}...")
        z_norm_df.write_parquet(zscore_path, compression="zstd", compression_level=22)

        size_mb = os.path.getsize(zscore_path) / (1024 * 1024)
        print(f"  Shape: {z_norm_df.shape}")
        print(f"  Output size: {size_mb:.2f} MB")

    print("\n" + "="*60)
    print("Done!")
    print("="*60)
    print("\nOutput files:")
    print(f"  - {importance_path.name}: Feature importance per variant")
    if not args.skip_features:
        print(f"  - {output_path.name}: Minimal single-cell features")
        print(f"  - {zscore_path.name}: Z-score normalized variant features")


if __name__ == "__main__":
    main()
