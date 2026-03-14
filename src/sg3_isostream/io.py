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
