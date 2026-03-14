# SG3_Isostream
### **Next-Generation Single-Cell Long-Read RNA-seq Analysis Suite**

**Developed by:** **Dr. Saminathan Sivaprakasham Murugesan** *Principal Investigator, Sumathi Genomics Suite* **Sumathi Genomics Co Ltd (Thailand)** | **Saminathan Industries Pte Ltd (Singapore)** **Support:** [lab@sumathigenomics.com](mailto:lab@sumathigenomics.com)

---

## 1. Introduction
The landscape of transcriptomics is shifting from gene-level quantification to isoform-level resolution. Traditional short-read RNAseq excels at counting total gene expression but often fails to resolve complex alternative splicing events. SingleCell Long-Read (SCLR) sequencing powered by Oxford Nanopore (ONT) and PacBio—offers the potential to sequence full-length transcripts.

However, SCLR data introduces unique computational challenges, primarily **high indel rates** and **significant transcript truncation**. **SG3_Isostream** is a Python native suite designed to rescue these truncated reads using a "Fuzzy-Boundary" recovery algorithm, enabling accurate isoform-level discovery within the `Scanpy` ecosystem.

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

import pysam
from collections import defaultdict
import tqdm

def parse_gtf(gtf_path):
    """
    Parses a GTF file to create a dictionary of gene models.
    Returns: { gene_id: { isoform_id: [(start, end), (start, end)] } }
    """
    gene_models = defaultdict(lambda: defaultdict(list))
    
    print(f"Parsing GTF: {gtf_path}...")
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if fields[2] != 'exon': continue
            
            # Extract attributes
            attrs = {x.split(' ')[0]: x.split(' ')[1].replace('"', '') 
                     for x in fields[8].split('; ') if x}
            
            gene_id = attrs.get('gene_id')
            transcript_id = attrs.get('transcript_id')
            
            # Append exon coordinates (0-based for compatibility with pysam)
            start, end = int(fields[3]) - 1, int(fields[4])
            gene_models[gene_id][transcript_id].append((start, end))
            
    # Sort exons for each transcript by start position
    for g_id in gene_models:
        for t_id in gene_models[g_id]:
            gene_models[g_id][t_id].sort()
            
    return gene_models

def load_bam_to_counts(bam_path, gene_models, assigner, cb_tag="CB"):
    """
    Streams a BAM file and assigns reads to isoforms.
    Returns: A nested dictionary of {cell_barcode: {isoform_id: count}}
    """
    counts = defaultdict(lambda: defaultdict(int))
    
    # Create an interval map for quick gene lookup
    # In production, use an IntervalTree for O(log n) lookup
    print(f"Processing BAM: {bam_path}...")
    
    with pysam.AlignmentFile(bam_path, "rb") as sam:
        for read in tqdm.tqdm(sam.fetch(), desc="Assigning Reads"):
            if read.is_unmapped or not read.has_tag(cb_tag):
                continue
                
            cb = read.get_tag(cb_tag)
            read_exons = read.get_blocks() # List of (start, end)
            
            # Simple overlap check: find which gene the read overlaps
            # This logic assumes the BAM is sorted/indexed
            ref_name = read.reference_name
            # (Note: In a full release, we match chromosome names between GTF and BAM)
            
            # Assignment logic
            assigned_iso = None
            for gene_id, isoforms in gene_models.items():
                # Preliminary check: does the read overlap the gene boundaries?
                gene_min = min(ex[0] for tx in isoforms.values() for ex in tx)
                gene_max = max(ex[1] for tx in isoforms.values() for ex in tx)
                
                if read.reference_start > gene_max or read.reference_end < gene_min:
                    continue
                
                assigned_iso = assigner.assign_read(read_exons, isoforms)
                if assigned_iso != "Unassigned_Fragment":
                    counts[cb][assigned_iso] += 1
                    break 
                    
    return counts
