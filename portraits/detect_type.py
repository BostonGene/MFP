import pandas as pd
from portraits.clustering import clustering_profile_metrics, clustering_profile_metrics_plot
from portraits.utils import read_gene_sets, ssgsea_formula, median_scale


def detect_type(data, threshold, scores):
    ser = data.loc[threshold].perc # here threshold and ser were added to the original code 
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

