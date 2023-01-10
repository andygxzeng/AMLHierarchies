## A cellular hierarchy framework for understanding heterogeneity and predicting drug response in AML
Publication: [Zeng et al, Nature Medicine 2022](https://www.nature.com/articles/s41591-022-01819-x)

This repository contains analysis notebooks and scripts corresponding to the main figures, as well as scRNA-seq data and deconvolution results used in the paper.  


**AML represents a caricature of normal hematopoietic development**, and this developmental process is distorted in different ways for different patients. Our study aimed to understand how leukemia cell hierarchies vary from patient to patient and how this relates to the functional, genomic, and clinical properties of each patient's disease.  

This analysis started with a focused [re-analysis of primitive AML cells](https://github.com/andygxzeng/AMLHierarchies/blob/main/Fig_1.1_LSPC_Reclustering_Characterization.ipynb) at the apex of leukemia cell hierarchies and applied deconvolution to understand how each of these primitive cell types relate to [functional LSC activity](https://github.com/andygxzeng/AMLHierarchies/blob/main/Fig_1.2_LSC_Analysis.ipynb). By profiling the leukemic hierarchy compositions of [over 1000 AML patients](https://github.com/andygxzeng/AMLHierarchies/blob/main/Fig_2.0_Cohort_Cluster_Survival.ipynb), we found that hierarchy composition was associated with [survival outcomes](https://github.com/andygxzeng/AMLHierarchies/blob/main/Fig_2.2_Hierarchy_Survival.ipynb), [genomic alterations](https://github.com/andygxzeng/AMLHierarchies/blob/main/Fig_2.1_Hierarchy_Genomic_Correlates.ipynb), and [disease relapse](https://github.com/andygxzeng/AMLHierarchies/blob/main/Fig_3.0_Relapse_Deconvolution.ipynb). 

Applied to drug screening data, we found that cells residing at different levels of the AML hierarchy differed in their [drug sensitivity profiles](https://github.com/andygxzeng/AMLHierarchies/blob/main/Fig_4.0_CellType_Drug_Correlations.ipynb) and trained simple gene expression scores to approximate hierarchy composition and [predict drug response](https://github.com/andygxzeng/AMLHierarchies/blob/main/Fig_4.1_LinClass7_PC2_Regression.ipynb). To apply this framework to drug development, we re-analyzed published pre-clinical studies from the literature to show how each drug treatment condition [affected cell type composition](https://github.com/andygxzeng/AMLHierarchies/blob/main/Fig_5.1a_Literature_Deconv_Composition_UMAP.ipynb) in these studies. Last, we showed that stratifiying patient samples based on hierarchy can robustly distinguish drug responders from non-responders in [patient-derived xenograft models](https://github.com/andygxzeng/AMLHierarchies/blob/main/Fig_6.0_Hierarchy_PDX_Response.ipynb). Together, this establishes a new framework for understanding AML heterogeneity with important implications for precision medicine efforts in AML. 



Deconvolution results are included in the "Data" directory, according to each analysis section. Due to large file sizes, re-annotated single cell RNA-seq data from AML patients (from van Galen et al) are hosted on AWS:   

**Re-annotated scAML data (from van Galen et al)**   
  - [Leukemia Stem and Progenitor Cells (LSPCs) only - scran normalized](https://amlhierarchies.s3.us-east-2.amazonaws.com/scRNA_analysis/scAML_LSPC_scran_reclustered.h5ad)  
  - [Malignant and immune cells - raw counts](https://amlhierarchies.s3.us-east-2.amazonaws.com/scRNA_analysis/scAML_Full_reannotated_rawcounts.h5ad)

**scAML TF regulon analysis (pySCENIC, malignant cells only)**
  - [TF Regulon Enrichment Scores](https://amlhierarchies.s3.us-east-2.amazonaws.com/scRNA_analysis/scAML_pySCENIC_regulon_scores_AUCell.csv)
  - [TF Regulon Information](https://amlhierarchies.s3.us-east-2.amazonaws.com/scRNA_analysis/scAML_pySCENIC_regulon_info.csv)
  
   
   
## AML Deconvolution Instructions
Through benchmarking experiments in our paper we have identified CIBERSORTx to perform best in deconvoluting AML data with our reference cell types. We have prepared two signature matrices for use in performing CIBERSORTx deconvolution on your TPM-normalized RNA-seq data. 

**CIBERSORTx Deconvolution (Malignant + Immune)**  
This signature matrix is comprised of 7 malignant cell types and 7 immune cell types and can be applied to any unsorted AML sample with infiltrating immune cells – we provide RNA-seq data from the TCGA cohort as an example dataset. 
  - [Full Signature Matrix](https://amlhierarchies.s3.us-east-2.amazonaws.com/Deconvolution/CIBERSORTx_scAML_Full_SignatureMatrix.txt)
  - [Full Single Cell Reference Sample](https://amlhierarchies.s3.us-east-2.amazonaws.com/Deconvolution/CIBERSORTx_scAML_Full_SingleCellReferenceSamp.txt)
  - [Example Dataset: TCGA Cohort](https://amlhierarchies.s3.us-east-2.amazonaws.com/Deconvolution/TCGA_LAML_TPM.txt)

**CIBERSORTx Deconvolution (Malignant only)**  
This signature matrix is comprised only of the 7 malignant cell types, with no immune populations. This can be applied to sorted AML samples or in experimental settings (e.g. cell lines, cultured primary samples, PDX models) – we provide data from sorted LSC fractions as an example dataset. 
  - [Malignant Signature Matrix](https://amlhierarchies.s3.us-east-2.amazonaws.com/Deconvolution/CIBERSORTx_scAML_Malignant_SignatureMatrix.txt)
  - [Malignant Single Cell Reference Sample](https://amlhierarchies.s3.us-east-2.amazonaws.com/Deconvolution/CIBERSORTx_scAML_Malignant_SingleCellReferenceSamp.txt)
  - [Example Dataset: Sorted AML fractions annotated for LSC activity](https://amlhierarchies.s3.us-east-2.amazonaws.com/Deconvolution/AML_LSC_fractions_TPM.txt)

To run CIBERSORTx, we recommend using the web portal (due to discrepancies in batch correction behaviour between the web portal and docker version) and applying deconvolution in Absolute mode using the provided Signature Matrix and Mixture (bulk) dataset, while applying S-mode batch correction using the provided single cell reference sample. Permutations are optional. 

After deconvolution, we recommend normalizing the malignant cell populations to 1 and projecting your samples onto the reference cohort distribution (TCGA, BEAT, Leucegene) using built-in functions from scanpy or Seurat. A simpler way to project your samples (if you have a small sample size) is to concatenate them with the reference samples, apply ComBat batch correction, and re-run PCA altogether. You can refer to our notebooks, particularly the [Relapse Deconvolution](https://github.com/andygxzeng/AMLHierarchies/blob/main/Fig_3.0_Relapse_Deconvolution.ipynb) notebook, for examples on projecting and analyzing new deconvolution data. 



**Signature enrichment analysis with AML cell type-specific genesets**  

If you prefer to perform signature scoring at the bulk level through GSVA or ssGSEA, we provide a gmt file with genesets specific to each AML cell type within the scRNA-seq data. Genesets for each AML cell type (LSPC-Quiescent, LSPC-Primed, LSPC-Cycling, GMP-like, ProMono-like, Mono-like, cDC-like) were generated by differential gene expression analysis through MAST, comparing each individual leukemic population against all other populations. When more than 250 DE genes were identified, genesets are restricted to the top 100 and top 250 DE genes for each population. Additional LSC genesets from Ng et al 2017 (pan-AML, identified through sorting and xenotransplantation) and Sommervaile et al 2009 (MLL-specific) are also provided.

The genesets can be found in the following directory [Data/AMLCellType_Genesets.gmt](https://raw.githubusercontent.com/andygxzeng/AMLHierarchies/main/Data/AMLCellType_Genesets.gmt)
