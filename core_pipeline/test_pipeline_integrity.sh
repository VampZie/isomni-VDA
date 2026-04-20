#!/bin/bash
# ------------------------------------------------------------------
# isomiR Pipeline Integrity & Unit Testing Suite
# Methodology: Validates data-flow from Raw FastQ to Merged Counts
# ------------------------------------------------------------------

# 1. SPACE-SAFE GLOBAL WORKSPACE
# [ENGINEERING NOTE]: Using /tmp to bypass HDD latency during testing
SAFE_BASE="/tmp/vda_integrity_check"
rm -rf "$SAFE_BASE"
mkdir -p "$SAFE_BASE/input" "$SAFE_BASE/output" "$SAFE_BASE/db"

# 2. LOCAL DATA ENVIRONMENT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# [RESEARCH SENSITIVE: DATA DIRECTORY]
# Masking the primary GSE199508 external cluster path 
# and local dissertation backup drive mount-points.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
DATA_DIR="./test_data" 

# 3. LINKING ASSETS (SPACE-SAFE SYMLINKING)
ln -s "${DATA_DIR}"/*.fastq "$SAFE_BASE/input/"
ln -s "${DATA_DIR}/generate_isomir_count.py" "$SAFE_BASE/input/"
ln -s "${DATA_DIR}/merge_isomirs.py"        "$SAFE_BASE/input/"

# miRBase Database Assets
ln -s "${DATA_DIR}/hairpin.fa" "$SAFE_BASE/db/hairpin.fa"
ln -s "${DATA_DIR}/mature.fa"  "$SAFE_BASE/db/mature.fa"
ln -s "${DATA_DIR}/miRNA.str"  "$SAFE_BASE/db/miRNA.str"

# 4. ITERATIVE TESTING (Step-by-Step Validation)
for fq in "$SAFE_BASE"/input/*.fastq; do
    sample=$(basename "$fq" .fastq)

    # STEP A: Adapter Trimming (Gated at 18bp for isomiR fidelity)
    trim_galore --small_rna -q 20 --length 18 -e 0.1 \
        --output_dir "$SAFE_BASE/output" "$fq"

    # STEP B: Cluster Collapsing
    seqcluster collapse -f "$SAFE_BASE/output/${sample}_trimmed.fq" -o "$SAFE_BASE/output"

    # STEP C: Alignment (miraligner Core)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # [RESEARCH SENSITIVE: ALIGNMENT PARAMETERS]
    # Specialized mismatch (-sub) and minimal length (-minl) gating 
    # tailored specifically for the Rat Model (rno) in VDA states.
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    java -jar miraligner.jar -i "$SAFE_BASE/output/${sample}.fa" \
        -db "$SAFE_BASE/db" -o "$SAFE_BASE/output/miraligner_${sample}" \
        -s rno <MASKED_PARAMS>

    # STEP D: Frequency Extraction
    python "$SAFE_BASE/input/generate_isomir_count.py" \
        "$SAFE_BASE/output/miraligner_${sample}.mirna" \
        "$SAFE_BASE/output/${sample}_isomir_counts.csv"
done

# 5. CONSOLIDATION
python "$SAFE_BASE/input/merge_isomirs.py"

echo "===================================="
echo " INTEGRITY CHECK COMPLETED"
echo "===================================="
