import pandas as pd
import argparse
import logging
import sys

# Configure professional logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def extract_isomirs(input_mirna, output_csv):
    """
    Refined logic for extracting isomiR frequencies from miraligner output.
    Captures 5' end shifts (t5) which are critical for non-canonical targeting 
    in Neurodegenerative disease models.
    """
    try:
        logging.info(f"Processing miraligner data: {input_mirna}")
        df = pd.read_csv(input_mirna, sep='\t')
        
        # Clean sequences: Remove added nucleotides generated during alignment
        df['seq_clean'] = df.apply(lambda r: r['seq'][:-len(r['add'])] if str(r['seq']).endswith(str(r['add'])) else r['seq'], axis=1)
        
        # Methodology: Aggregating by 5' shift (t5) to find novel variants
        # Lowercase or 0 tags in t5 indicate variations from canonical start sites
        grouped = df.groupby(['mir', 'precursor', 't5'], as_index=False).agg({
            'freq': 'sum',
            'seq_clean': 'first'
        })
        
        # Rename MIR to include shift metadata for downstream dMiso integration
        grouped['mir_variant'] = grouped.apply(lambda r: f"{r['t5']}_{r['mir']}" if str(r['t5']) != '0' else str(r['mir']), axis=1)
        
        # Final formatting
        final_df = grouped[['mir_variant', 'freq', 'precursor', 'seq_clean']]
        final_df.to_csv(output_csv, index=False)
        logging.info(f"Successfully exported {len(final_df)} isomiR variants to {output_csv}")

    except Exception as e:
        logging.error(f"Failed to process isomiR data: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='VDA isomiR Extraction Tool')
    parser.add_argument('--input', required=True, help='Input miraligner .mirna file')
    parser.add_argument('--output', required=True, help='Output path for count CSV')
    args = parser.parse_args()
    extract_isomirs(args.input, args.output)
