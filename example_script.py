# Python3.7
import pandas as pd

from portraits.clustering import clustering_profile_metrics, clustering_profile_metrics_plot
from portraits.utils import read_gene_sets, ssgsea_formula, median_scale

# Example script

# Read signatures
gmt = read_gene_sets('signatures.gmt')  # GMT format like in MSIGdb

# Read expressions
exp = pd.read_csv('expression.tsv', sep='\t', index_col=0)  # log2+1 transformed; Genes in columns

# Calc signature scores
signature_scores = ssgsea_formula(exp, gmt)

# Scale signatures
signature_scores_scaled = median_scale(signature_scores)

# Check the clustering within a range of 30 to 65% similarity.
# >65% - usually graph is not connected; <30% - unreasonable correlation
clustering_metrics = clustering_profile_metrics(signature_scores_scaled, threshold_mm=(.3, .65), step=.01)

# Visualize the partitions
clustering_profile_metrics_plot(clustering_metrics)

# Then select the best threshold using one ore more metrics.
best_threshold = '0.51'


def detect_type(ser, scores):
    cmeans = pd.DataFrame({cg: scores.loc[samps.index].mean() for cg, samps in ser.groupby(ser)})
    mapper = {}
    deltas = (cmeans.loc[['Angiogenesis', 'Endothelium', 'CAF', 'Matrix', 'Matrix_remodeling']].mean() -
              cmeans.loc[['MHCII', 'Antitumor_cytokines', 'Coactivation_molecules',
                          'B_cells', 'NK_cells', 'Checkpoint_inhibition',
                          'Effector_cells', 'T_cells', 'Th1_signature',
                          'T_cell_traffic', 'MHCI']].mean()).sort_values()

    mapper[deltas.index[-1]] = 'F'  # That's fibrotic
    mapper[deltas.index[0]] = 'IE'  # Immune enriched, non-fibrotic
    cmeans.pop(deltas.index[-1])
    cmeans.pop(deltas.index[0])

    deltas = (cmeans.loc[['Angiogenesis', 'Endothelium', 'CAF', 'Matrix', 'Matrix_remodeling',
                          'Protumor_cytokines', 'Neutrophil_signature', 'Granulocyte_traffic',
                          'Macrophages', 'Macrophage_DC_traffic', 'MDSC_traffic', 'MDSC',
                          'Th2_signature', 'T_reg_traffic', 'Treg', 'M1_signatures', 'MHCII',
                          'Antitumor_cytokines', 'Coactivation_molecules', 'B_cells', 'NK_cells',
                          'Checkpoint_inhibition', 'Effector_cells', 'T_cells', 'Th1_signature',
                          'T_cell_traffic', 'MHCI', 'EMT_signature']].mean() -
              cmeans.loc['Proliferation_rate']).sort_values()

    mapper[deltas.index[-1]] = 'IE/F'  # Immune enriched & fibrotic
    mapper[deltas.index[0]] = 'D'  # Desert
    return ser.map(mapper).rename('MFP')


# Detect cluster types
final_clusters = detect_type(clustering_metrics[best_threshold], signature_scores_scaled)

# Output the clusters
final_clusters.to_csv('final_clusters.tsv', sep='\t', index=True)
