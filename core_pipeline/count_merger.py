import pandas as pd
import glob
import argparse
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def merge_counts(input_dir, output_file):
    """
    Groups individual sample isomiR counts into a single expression matrix.
    Uses Inner Join to ensure consistent reporting across the cohort.
    """
    files = glob.glob(os.path.join(input_dir, "*_isomir_counts.csv"))
    if not files:
        logging.warning("No isomiR count files found in input directory.")
        return

    merged_df = None
    logging.info(f"Merging {len(files)} samples...")

    for file in files:
        sample_id = os.path.basename(file).replace("_isomir_counts.csv", "")
        df = pd.read_csv(file)
        
        # Handle potential duplicates in MIR variant IDs
        if df['mir_variant'].duplicated().any():
            df['mir_variant'] = df.apply(
                lambda r: f"{r['mir_variant']}_{r['precursor']}" if df['mir_variant'].duplicated(keep=False)[r.name] else r['mir_variant'],
                axis=1
            )
            
        df = df[['mir_variant', 'freq']].rename(columns={'freq': sample_id})
        
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on='mir_variant', how='inner')

    if merged_df is not None:
        merged_df.to_csv(output_file, index=False)
        logging.info(f"Unified expression matrix saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='IsomiR Count Matrix Generator')
    parser.add_argument('--dir', default='.', help='Input directory containing sample CSVs')
    parser.add_argument('--output', default='merged_isomir_counts.csv', help='Output matrix file name')
    args = parser.parse_args()
    merge_counts(args.dir, args.output)
