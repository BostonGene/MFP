# Conserved pan-cancer microenvironment subtypes.

## Introduction
Tumor microenvironment (TME) plays a significant role in clinical outcomes and response to antineoplastic therapy. By exerting pro-tumorigenic and anti-tumorigenic actions, tumor-infiltrating immune cells can profoundly influence tumor progression and affect the success of anti-cancer treatments. Cancer-associated fibroblasts (CAFs), as well as angiogenic signals from stromal cells, have also been shown to affect outcomes in cancer patients. 


BostonGene has compiled a curated list of 29 functional genes and uses single-sample Gene Set Enrichment Analysis (ssGSEA) scores of their expression levels to classify a tumor sample into one of the TME subtypes (Fig. 1). Materials provided in this repository will help to identify the TME type of an input sample.

![image2020-11-4_18-20-47](https://user-images.githubusercontent.com/127855909/228009303-964b1147-0f42-4361-819b-bc22be9ccd97.png)

<p align="center">Figure 1. Types of tumor microenvironment.</p>

<p align="center">A brief description of each type of TME and its graphical interpretation can be found in Fig. 2.</p>

![image2020-11-4_18-12-26](https://user-images.githubusercontent.com/127855909/228009221-3fe09cc9-a30a-4d3f-aa4b-3641c6278f7e.png)

<p align="center">Figure 2. TME types brief description.</p>

## Citation
If software, data, and/or website are used in your publication, please cite [Bagaev A et al. Conserved pan-cancer microenvironment subtypes predict response to immunotherapy. Cancer Cell. 2021 Jun 14;39(6):845-865](https://www.cell.com/cancer-cell/fulltext/S1535-6108(21)00222-1#articleInformation) and make a reference to this repository.


For more information visit [BostonGene’s Scientific portal](https://science.bostongene.com/tumor-portrait/)


## Setup
Set up your environment according to the requirements in the description of the Setup.md file in the repository.


If your environment is already set up accordingly, clone the Github repository to start your analysis.



    git clone https://github.com/BostonGene/MFP.git


## Implementation overview
***Note: The example analysis is performed for a cohort. Do not perform TME classification analyses for just one sample.***


The analysis workflow is presented in the diagram below, highlighting the main steps and logical elements of the notebook.

![TME  workflow](https://github.com/BenjaminSargsyan/Preview_of_TME/assets/127855909/037999b3-4f43-4770-bf5a-2e7f5ca9b909)


The analysis comprises two notebooks: [TME_Classification.ipynb](TME_Classification.ipynb) and [GEO_data_retrieval.ipynb](GEO_data_retrieval.ipynb).

TME_Classification.ipynb serves as the primary component of the analysis, encompassing various sections for data processing and classification. It includes an additional section allowing users to incorporate their own reference cohort for classification purposes.

For users intending to download their data from GEO before initiating the analyses, the GEO_retrieval.ipynb notebook is available. This notebook provides a dedicated framework for data retrieval specifically from GEO, enabling seamless integration into the subsequent analysis pipeline.


The pipeline consists of several nodes that correspond to each other where some of them are optional(depends on user choice):

* Data preparation
  * Data retrieval (and normalization) 
* Quality Check (QC)
  * Batch detection
  * Outliers detection
  * Data distribution check for data quality
* Classification
  * Reading the expressions.tsv file to get the gene expression matrix of the samples of interest
  * Getting the reference gene signatures expressions matrix (TCGA cohort is set by default and can be changed to the path to your reference cohort)
  * Identifying the TME subtype of the sample/samples of interest by comparing their ssGSEA score to the ssGSEA scores of the reference cohort
  * Giving an output .tsv file with the sample/samples TME subtype
* Clusterization of a reference cohort (optional; we recommend using the default TCGA cohort)
  * Getting a reference cohort input
  * Identifying the TME subtype for each sample based on its ssGSEA score
  * Getting an output .tsv table with the sample subtypes in the reference cohort

**Note: It is recommended to use the default TCGA cohort to avoid possible problems during analysis.**


You can also access the [Data_processing_methods.ipynb](Data_processing_methods.ipynb) notebook, which provides additional information on how the batch correction and outlier detection analyses were done in the article.
