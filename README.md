HEAD
# VDA-miRNA-Collaborative 🧬
### *A High-Resolution Framework for isomiR-mRNA Dynamics in Vascular Dementia*

---

## 📋 Overview
This repository contains the multi-omics computational framework developed for the comparative analysis of novel and known miRNAs in **Vascular Dementia (VDA)**. Unlike traditional pipelines, this framework explicitly accounts for **isomiR (miRNA isoform)** dynamics—capturing 5' and 3' end variations that significantly shift regulatory targets in neurodegenerative states.

**This project is a collaborative research effort** focusing on identifying non-canonical miRNA-mRNA interaction networks to uncover novel biomarkers for VDA progression.

---

## 🛠️ Tech Stack & Dependencies
- **Orchestration:** Shell / Nextflow (Ready)
- **Languages:** Python (Pandas, BioPython), R (DESeq2, clusterProfiler, ReactomePA)
- **Tools:** miraligner, miRanda, dMiso (Isomer-aware prediction)
- **Environment:** Docker / Conda

---

## 🚀 Pipeline Architecture

### Phase 1: Pre-processing & isomiR Discovery
- **QC & Trimming:** Small-RNA specific adapter removal using `trim_galore`.
- **Mapping:** Mapping reads to miRBase (v22) using `miraligner` to extract non-canonical shifts.
- **Aggregation:** Custom frequency aggregation logic to identify dominant isomiR species.

### Phase 2: Differential Expression (DE)
- **mRNA/miRNA Analysis:** Concurrent DE analysis using `DESeq2` with `apeglm` LFC shrinkage to handle low-frequency variant noise.

### Phase 3: Regulatory Network Validation
- **Target Prediction:** Cross-referencing miRanda (energetic stability) with `dMiso` (Isomer-specific thermodynamic modeling).
- **Network Construction:** Correlation of miRNA-isomiR abundance with target mRNA suppression.

### Phase 4: Functional Pathfinding
- **Enrichment:** Automated GO, KEGG, and **Reactome** pathway analysis to identify neuro-inflammatory signals.

---

## 📂 Repository Structure
```text
├── core_pipeline/           # isomiR extraction and frequency merging
├── differential_expression/ # DESeq2 workflows for mRNA and miRNA
├── interaction_prediction/  # Sequence pair generation and dMiso orchestration
├── functional_analysis/     # Pathway enrichment and visualization (R)
└── docs/                    # Methodology and environment specifications
```

---

## 📈 Key Outcomes
- Identification of specific **5' shifts** in neuroprotective miRNAs that alter seed-region complementary.
- Discovery of novel VDA-specific miRNA isoforms through rigorous frequency filtering.
- Validated regulatory flips in synaptic plasticity genes through thermodynamic modeling.

---

## 👥 Authors & Collaboration
This project represents a collaborative effort in the field of Precision Neurology. 
- **Lead Developer:** Vidit Zainith ([@VampZie](https://github.com/VampZie))
- **Status:** Manuscript in Communication / Dissertation Objective 1.

---
*Developed as part of the M.Sc. Bioinformatics Dissertation Pipeline.*
=======
# isomni-VDA
Novel Isoform Analysis
ab29c5a9de5cc9cdf723dc8fd82f132d6bb948a7
