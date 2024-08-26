# scRNA-seq_Tutorial

This repository contains the pre-processing and analysis steps for scRNA-seq (single cell RNA-seq) data from one cell type or one experimental condition.

For this tutorial, we will be analyzing the a dataset of Non-Small Cell Lung Cancer Cells (NSCLC) freely available from 10X Genomics (https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_nsclc_5gex), using the Seurat R package (http://satijalab.org/seurat/), a popular and powerful set of tools to conduct scRNA-seq analysis in R. In this dataset, there are 7802 single cells that were sequenced on the Illumina NovaSeq 6000. Please note this tutorial borrows heavily from Seurat’s tutorials (https://satijalab.org/seurat/vignettes.html), so feel free to go through them in more detail.

The data for Non-Small Cell Lung Cancer Cells (NSCLC) is freely available from 10X Genomics (https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_nsclc_5gex).

The repository has scRNA.R (a R script for the analysis of scRNA data), which has been incorporated in scRNA.sh (which has pre-processing steps of the original data).

The Rscript has all the processes and steps defined with a "#", for the user's convenience. 

All the image files (heatmaps, plots, and others) generated from the R script have been provided here.

# Integration of scRNA-seq data from two different conditions

Integration of single-cell sequencing datasets, for example across experimental batches, donors, or conditions, is often an important step in scRNA-seq workflows.

The scRNA_integrated.R script performs the integration of scRNA-seq data from human PBMC from two different conditions.

The object "ifnb" contains data from human PBMC from two conditions, interferon-stimulated and control cells (stored in the stim column in the object metadata). We will aim to integrate the two conditions together, so that we can jointly identify cell subpopulations across datasets, and then explore how each group differs across conditions.

The Rscript has all the processes and steps defined with a "#", for the user's convenience. 

All the image files (heatmaps, plots, and others) generated from the R script have been provided here.
