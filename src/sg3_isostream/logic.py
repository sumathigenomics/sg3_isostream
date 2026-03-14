import numpy as np

class SG3_Assigner:
    def __init__(self, tolerance=50, junction_weight=0.8):
        """
        Engine for SCLR-seq read assignment.
        tolerance: bp slack for exon boundaries.
        junction_weight: Importance of internal splice sites vs transcript ends.
        """
        self.tolerance = tolerance
        self.j_weight = junction_weight
        self.e_weight = 1.0 - junction_weight

    def _get_junctions(self, exons):
        """Converts exon coordinates into unique splice-junction pairs."""
        return set([(exons[i][1], exons[i+1][0]) for i in range(len(exons)-1)])

    def recovery_score(self, read_exons, ref_exons):
        """
        Calculates the probability that a read belongs to a reference transcript.
        Focuses on internal junctions to rescue 5' or 3' truncated reads.
        """
        read_junc = self._get_junctions(read_exons)
        ref_junc = self._get_junctions(ref_exons)
        
        # 1. Junction Score: Do the splice sites match?
        if not read_junc:
            # Single-exon read logic: check if read is within ref boundaries
            j_score = 1.0 if (read_exons[0][0] >= ref_exons[0][0] - self.tolerance and 
                             read_exons[0][1] <= ref_exons[0][1] + self.tolerance) else 0.0
        else:
            matches = len(read_junc.intersection(ref_junc))
            j_score = matches / len(read_junc)
            
        # 2. Boundary Match: Is the read a subset of the reference?
        is_subset = (read_exons[0][0] >= ref_exons[0][0] - self.tolerance and 
                     read_exons[-1][1] <= ref_exons[-1][1] + self.tolerance)
        
        return (j_score * self.j_weight) + (float(is_subset) * self.e_weight)

    def assign_read(self, read_exons, gene_models):
        """Maps a noisy read to the best matching isoform in a gene model."""
        best_iso, best_score = None, 0
        for iso_id, ref_exons in gene_models.items():
            s = self.recovery_score(read_exons, ref_exons)
            if s > best_score:
                best_score, best_iso = s, iso_id
        
        return best_iso if best_score > 0.75 else "Unassigned_Fragment"
