# The repository contains code and materials linked with work "Integrated tumor exome and transcriptome analyses reveal conserved pan-cancer microenvironment subtypes predictive of response to immunotherapy"

### Navigation

1. Renormalization form CEL files is shown in <a href="https://nbviewer.jupyter.org/github/BostonGene/MFP/blob/master/From_cell_files.ipynb">From_cell_files.ipynb</a>
1. Signatures in gmt format could be found in signatures/ folder 
1. Outliers detection and batch detection & correction are shown in <a href="https://nbviewer.jupyter.org/github/BostonGene/MFP/blob/master/Methods_Description_-_Batch_correction.ipynb">Methods_Description_-_Batch_correction.ipynb</a> (Melanoma data)
1. Clusters generation from scaled signatures is shown in <a href="clustering_example.py">clustering_example.py</a>. The process is the same for melanoma and pancan analysis
1. Finally, having a dataset with scaled signatures and known clusters we can classify another datasets using <a href="classification_example.py">classification_example.py</a>


_.ipynb files could by opened at https://nbviewer.jupyter.org/ or downloaded as HTML files from upstream_html folder_


# MFP - Molecular Functional Portrait

The Molecular Functional (MF) Portrait is a planetary schematic representation of integrated molecular and functional characteristics of a tumor and its microenvironment based on whole exome (WES) and RNAseq analysis. The MF Portrait depicts the prevalence of malignant and tumor microenvironment (TME) cell populations, and the activity of tumor promoting and suppressing processes. The portrait includes qualitative and quantitative descriptions of modules built based on the expression of BostonGene curated 29 gene expression signatures (reference manuscript and table), with the size of each module corresponding to the intensity of the normalized single-sample gene set enrichment analysis (GSEA) score, and the colors denoting pro- (red) or anti-cancer activity (blue). In the MF Portrait, potentially targetable genomic alterations and TME processes such as cell proliferation, EMT, angiogenesis, and anti-tumor immunity are listed. The relevant TME subtypes (reference our paper) for each tumor, identified based on BostonGene classification platform, is are integrated into the MF Tumor Portrait to aid in the rational design of therapeutic strategies.


![Graphical abstract](img/Abstract.svg?raw=true "Molecular Functional Portrait")

## TME subtypes include:

* Immune-enriched, fibrotic (IE/F)
* Immune-enriched (IE)
* Fibrotic (F)
* Depleted (D)

![mfp_characteristics](img/mfp_characteristics.png?raw=true "MFP characteristics")

Visual tool available at https://science.bostongene.com/tumor-portrait/<br>

© 2020 BostonGene Corporation.