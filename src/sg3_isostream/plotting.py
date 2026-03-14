import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def plot_sashimi(isocell, gene_id, cluster_key, clusters_to_plot):
    """
    Plots read coverage and splice-junction arcs for specific cell clusters.
    """
    fig, axes = plt.subplots(len(clusters_to_plot), 1, 
                             figsize=(12, 5 * len(clusters_to_plot)), 
                             sharex=True)
    
    if len(clusters_to_plot) == 1: axes = [axes]

    # In a real run, these coordinates would come from the GTF/isocell metadata
    # Mocking coordinates for the gene
    gene_start, gene_end = 5000, 8000
    
    for i, cluster_name in enumerate(clusters_to_plot):
        ax = axes[i]
        
        # 1. Draw the 'Coverage Mountain' (mock data for illustration)
        x_coords = np.linspace(gene_start, gene_end, 1000)
        # Simulate read depth
        y_depth = np.random.normal(10, 2, 1000) + (np.sin(x_coords/100) * 5)
        y_depth = np.clip(y_depth, 0, None)
        
        ax.fill_between(x_coords, y_depth, color='#3498db', alpha=0.3)
        ax.plot(x_coords, y_depth, color='#2980b9', lw=1)

        # 2. Draw Junction Arcs (Splice Evidence)
        # Mocking two junctions: Exon 1-2 and Exon 2-3
        junctions = [(5500, 6200, 120), (6500, 7500, 85)] 
        
        for start, end, count in junctions:
            center = (start + end) / 2
            width = end - start
            height = count / 5 # Scale height by read count
            
            # Create a parabolic arc using PathPatch
            verts = [(start, 0), (center, height), (end, 0)]
            codes = [patches.Path.MOVETO, patches.Path.CURVE3, patches.Path.CURVE3]
            path = patches.Path(verts, codes)
            patch = patches.PathPatch(path, edgecolor='#e74c3c', facecolor='none', 
                                     lw=np.log1p(count), alpha=0.7)
            
            ax.add_patch(patch)
            ax.text(center, height + 0.5, f"n={count}", ha='center', fontsize=9, color='#c0392b')

        ax.set_title(f"Cluster: {cluster_name}", loc='left', fontweight='bold')
        ax.set_ylabel("Read Depth")
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.xlabel(f"Genomic Coordinates on {gene_id}")
    plt.tight_layout()
    return fig
