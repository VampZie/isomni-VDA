import pandas as pd
from Bio import SeqIO
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def generate_pairs(mirna_fasta, utr_fasta, interaction_csv, output_tsv):
    """
    Constructs miRNA-UTR sequence pairs for dMiso prediction.
    Integrates expression-based interaction hits with high-resolution 
    sequence database metadata.
    """
    logging.info("Building sequence lookup maps...")
    
    # Load miRNA sequences
    mirna_map = {rec.id: str(rec.seq) for rec in SeqIO.parse(mirna_fasta, "fasta")}
    
    # Load UTR sequences (keyed by Gene Symbol)
    utr_map = {}
    for rec in SeqIO.parse(utr_fasta, "fasta"):
        gene_id = rec.id.split("::")[0]
        utr_map.setdefault(gene_id, []).append(str(rec.seq))

    # Load predicted interactions (e.g. from miRanda)
    interactions = pd.read_csv(interaction_csv)
    interactions.columns = ["miRNA", "Gene", "Score", "Energy"]

    rows = []
    missing_mir = 0
    missing_utr = 0

    logging.info(f"Processing {len(interactions)} interaction candidates...")
    
    for _, row in interactions.iterrows():
        mir_id = row["miRNA"]
        gene_id = row["Gene"]

        if mir_id not in mirna_map:
            missing_mir += 1
            continue
        if gene_id not in utr_map:
            missing_utr += 1
            continue

        mir_seq = mirna_map[mir_id]
        for utr_seq in utr_map[gene_id]:
            rows.append([mir_id, gene_id, mir_seq, utr_seq])

    # Save to TSV (dMiso compatible)
    out_df = pd.DataFrame(rows, columns=["miRNA_ID", "Gene_ID", "miRNA_sequence", "UTR_sequence"])
    out_df.to_csv(output_tsv, sep="\t", index=False)
    
    logging.info(f"Pair generation complete. Total pairs: {len(out_df)}")
    logging.info(f"Skipped: {missing_mir} miRNAs missing sequence, {missing_utr} Genes missing UTR.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='dMiso Pair Input Generator')
    parser.add_argument('--mirna', required=True, help='miRNA FASTA file')
    parser.add_argument('--utr', required=True, help='UTR FASTA file')
    parser.add_argument('--hit', required=True, help='miRanda/Interaction CSV')
    parser.add_argument('--output', default='dmiso_pairs.tsv', help='Output TSV')
    args = parser.parse_args()
    generate_pairs(args.mirna, args.utr, args.hit, args.output)
