import pandas as pd
import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix

class IsoCell(AnnData):
    """
    The SG3_Isostream core data structure. 
    Maintains transcript-level sparsity while providing gene-level layers.
    """
    def __init__(self, transcript_matrix, obs, var, tx_to_gene_map):
        super().__init__(X=transcript_matrix, obs=obs, var=var)
        self.uns['tx_map'] = tx_to_gene_map
        self._collapse_to_genes()

    def _collapse_to_genes(self):
        """Sum isoforms to generate standard gene-level expression matrix."""
        genes = sorted(list(set(self.uns['tx_map'].values())))
        gene_to_idx = {g: i for i, g in enumerate(genes)}
        
        mapping_matrix = np.zeros((self.n_vars, len(genes)))
        for i, tx_id in enumerate(self.var_names):
            gene_id = self.uns['tx_map'].get(tx_id, "Unknown")
            mapping_matrix[i, gene_to_idx[gene_id]] = 1
            
        self.layers['gene_counts'] = self.X @ csr_matrix(mapping_matrix)
        self.uns['gene_names'] = genes
