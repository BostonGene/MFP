# Python3.7
import pandas as pd

from portraits.classification import KNeighborsClusterClassifier
from portraits.utils import read_gene_sets, ssgsea_formula, median_scale

# Example of how to classify another cohort having clusters on TCGA for example

# Load Training cohort with known MFP labels
TCGA_signature_scores_scaled = pd.read_csv('TCGA_signatures.tsv', sep='\t', index_col=0)  # Signatures in columns
TCGA_annotation = pd.read_csv('TCGA_annotation.tsv', sep='\t', index_col=0)  # Contains MFP cluster labels in MFP column
#  Fit the model
MODEL = KNeighborsClusterClassifier(norm=False, scale=False, clip=2, k=35).fit(TCGA_signature_scores_scaled,
                                                                               TCGA_annotation.MFP)

# Load the cohort of interest
# Read signatures
gmt = read_gene_sets('signatures.gmt')  # GMT format like in MSIGdb

# Read expressions
exp = pd.read_csv('expression.tsv', sep='\t', index_col=0)  # log2+1 transformed; Genes in columns

# Calc signature scores
signature_scores = ssgsea_formula(exp, gmt)

# Scale signatures
signature_scores_scaled = median_scale(signature_scores)

# Predict clusters
cluster_labels = MODEL.prtedict(signature_scores_scaled[MODEL.X.columns])

# Output the clusters
cluster_labels.to_csv('predicted_clusters.tsv', sep='\t', index=True)
