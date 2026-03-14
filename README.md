SG3_Isostream
Next-Generation Single-Cell Long-Read RNA-seq Analysis Suite
Developed by: Dr. Saminathan Sivaprakasham Murugesan Principal Investigator, Sumathi Genomics Suite Sumathi Genomics Co Ltd (Thailand) | Saminathan Industries Pte Ltd (Singapore) Support: lab@sumathigenomics.com

1. Introduction
The landscape of transcriptomics is shifting from gene-level quantification to isoform-level resolution. Traditional short-read RNA-seq excels at counting total gene expression but often fails to resolve complex alternative splicing events. Single-Cell Long-Read (SCLR) sequencing—powered by Oxford Nanopore (ONT) and PacBio—offers the potential to sequence full-length transcripts.

However, SCLR data introduces unique computational challenges, primarily high indel rates and significant transcript truncation. SG3_Isostream is a Python-native suite designed to rescue these truncated reads using a "Fuzzy-Boundary" recovery algorithm, enabling accurate isoform-level discovery within the Scanpy ecosystem.

2. Core Innovations
Weighted Junction-Matching (Fuzzy Recovery)
Standard tools often discard truncated reads. SG3_Isostream introduces a Recovery Score, prioritizing internal splice-junction consistency (weighted at 80%) over noisy transcript ends (20%). This allows the tool to "rescue" fragments and assign them to the most probable parent isoform.

Compositional Differential Isoform Usage (DIU)
Beyond simple fold-change, SG3_Isostream utilizes a statistical framework to detect Isoform Switching—biological events where total gene expression remains static, but the transcript ratio shifts significantly between cell clusters.

High-Fidelity Visualizations
SG3-Sashimi Plots: Detailed representations of read coverage and splice-junction arcs.

Iso-DotPlots: Multi-dimensional view of gene expression (size) and isoform ratio (color).

3. Installation
To install the development version from source:

Bash

git clone https://github.com/sumathigenomics/SG3_Isostream.git
cd SG3_Isostream
pip install .
Requirements:

Python >= 3.9

Anndata, Scanpy, Pysam, Matplotlib, Scipy, Pandas, Numpy

4. Quick Start Guide
SG3_Isostream is designed to be intuitive for users of the Scanpy framework.

Python

import sg3_isostream as sg3

# 1. Initialize data (Compatible with AnnData/Scanpy)
# counts: cell x transcript matrix
isodata = sg3.IsoCell(counts, obs, var, tx_to_gene_map)

# 2. Identify Significant Isoform Switching
# Finds genes where isoforms switch between 'Tumor' and 'Normal' clusters
switches = sg3.find_isoform_switches(
    isodata, 
    cluster_a_mask=isodata.obs['cell_type'] == 'Tumor',
    cluster_b_mask=isodata.obs['cell_type'] == 'Normal'
)

# 3. Generate Publication-Quality Sashimi Plot
sg3.plot_sashimi(isodata, gene_id="CD44", clusters_to_plot=["Tumor", "Normal"])
5. Directory Structure
Plaintext

SG3_Isostream/
├── src/
│   └── sg3_isostream/
│       ├── __init__.py     # API Exposure
│       ├── core.py         # IsoCell Class
│       ├── logic.py        # Fuzzy Assigner & Recovery Score
│       ├── stats.py        # DIU Math
│       ├── io.py           # BAM/GTF Data Loaders
│       └── plotting.py     # Sashimi & DotPlots
├── tests/                  # Unit tests and Mock data
├── pyproject.toml          # Build configuration
└── README.md               # Documentation
6. Citation
If you use this software in your research, please cite it as follows:

APA Style:

Sivaprakasham Murugesan, S. (2026). SG3_Isostream: A specialized Python suite for fuzzy-boundary recovery and isoform-level analysis in single-cell long-read RNA sequencing (Version 1.0.0) [Computer software]. Sumathi Genomics Co Ltd. https://github.com/sumathigenomics/SG3_Isostream

7. License
This project is licensed under the MIT License - see the LICENSE file for details.
