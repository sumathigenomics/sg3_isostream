# SG3_Isostream
### **Next-Generation Single-Cell Long-Read RNA-seq Analysis Suite**

**Developed by:** **Dr. Saminathan Sivaprakasham Murugesan** *Principal Investigator, Sumathi Genomics Suite* **Sumathi Genomics Co Ltd (Thailand)** | **Saminathan Industries Pte Ltd (Singapore)** **Support:** [lab@sumathigenomics.com](mailto:lab@sumathigenomics.com)

---

## 1. Introduction
The landscape of transcriptomics is shifting from gene-level quantification to isoform-level resolution. Traditional short-read RNAseq excels at counting total gene expression but often fails to resolve complex alternative splicing events. SingleCell Long-Read (SCLR) sequencing—powered by Oxford Nanopore (ONT) and PacBio—offers the potential to sequence full-length transcripts.

However, SCLR data introduces unique computational challenges, primarily **high indel rates** and **significant transcript truncation**. **SG3_Isostream** is a Python-native suite designed to rescue these truncated reads using a "Fuzzy-Boundary" recovery algorithm, enabling accurate isoform-level discovery within the `Scanpy` ecosystem.

## 2. Core Innovations

### **Weighted Junction-Matching (Fuzzy Recovery)**
Standard tools often discard truncated reads. SG3_Isostream introduces a **Recovery Score**, prioritizing internal splice-junction consistency (weighted at 80%) over noisy transcript ends (20%). This allows the tool to "rescue" fragments and assign them to the most probable parent isoform.

### **Compositional Differential Isoform Usage (DIU)**
Beyond simple fold-change, SG3_Isostream utilizes a statistical framework to detect **Isoform Switching**—biological events where total gene expression remains static, but the transcript ratio shifts significantly between cell clusters.

### **High-Fidelity Visualizations**
* **SG3-Sashimi Plots:** Detailed representations of read coverage and splice-junction arcs.
* **Iso-DotPlots:** Multi-dimensional view of gene expression (size) and isoform ratio (color).

---

## 3. Installation

To install the development version from source:

```bash
git clone [https://github.com/sumathigenomics/SG3_Isostream.git](https://github.com/sumathigenomics/SG3_Isostream.git)
cd SG3_Isostream
pip install .

import sg3_isostream as sg3

# 1. Initialize data (Compatible with AnnData/Scanpy)
isodata = sg3.IsoCell(counts, obs, var, tx_to_gene_map)

# 2. Identify Significant Isoform Switching
switches = sg3.find_isoform_switches(
    isodata, 
    cluster_a_mask=isodata.obs['cell_type'] == 'Tumor',
    cluster_b_mask=isodata.obs['cell_type'] == 'Normal'
)

# 3. Generate Publication-Quality Sashimi Plot
sg3.plot_sashimi(isodata, gene_id="CD44", clusters_to_plot=["Tumor", "Normal"])

