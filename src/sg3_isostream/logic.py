import numpy as np

class SG3_Assigner:
    def __init__(self, tolerance=50, junction_weight=0.8):
        self.tolerance = tolerance
        self.j_weight = junction_weight
        self.e_weight = 1.0 - junction_weight

    def get_junctions(self, exons):
        return set([(exons[i][1], exons[i+1][0]) for i in range(len(exons)-1)])

    def recovery_score(self, read_exons, ref_exons):
        """Calculates match probability between noisy read and reference."""
        read_junc = self.get_junctions(read_exons)
        ref_junc = self.get_junctions(ref_exons)
        
        # Junction Match (The strongest biological evidence)
        j_score = len(read_junc.intersection(ref_junc)) / len(read_junc) if read_junc else 0.0
        
        # Boundary Match (Accounting for truncation)
        is_contained = (read_exons[0][0] >= ref_exons[0][0] - self.tolerance and 
                        read_exons[-1][1] <= ref_exons[-1][1] + self.tolerance)
        
        return (j_score * self.j_weight) + (float(is_subset) * self.e_weight)

    def assign_read(self, read_exons, gene_models):
        best_iso, best_score = None, 0
        for iso_id, ref_exons in gene_models.items():
            s = self.recovery_score(read_exons, ref_exons)
            if s > best_score:
                best_score, best_iso = s, iso_id
        return best_iso if best_score > 0.75 else "Unassigned_Fragment"
