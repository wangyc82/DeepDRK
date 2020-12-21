# DeepDR
**D**eep learning **D**rug **R**esponses.

# Status

Active development

# Introduction

Accurately predicting the response of a cancer patient to a therapeutic agent is core goal of precision oncology. Among various obstacles hindering clinical translation, the lack of effective multimodal and multi-source data integration methods has become a bottleneck. DeepDR provides a systematic way to predict drug response from different resources, diverse cancer types, and various chemical compounds, based on kernel-based integration of multi-source data from genomics,
transcriptomics, epigenomics, and chemical properties of compounds.

# Usage

1. Installation

Before installing DeepDR please make sure you have installed R, and Rscript is available in your system path ($PATH).

git clone https://github.com/wangyc82/DeepDRv1

2. Preparing the input files

The input files includes the omics profile for cancer cells, including mutation at gene level, copy number at gene level, DNA methylation at gene level, and gene expression, and the drug properties, including chemical properties generated from chemical structures, and target proteins.

All these omics profiles must represent same set of cancer cells, that is, you have to prepare these omics data for same set of cancer cells, and let them is a R list data, such as cell_tst in example.data.

The chemical properties are the molecular fingerprints based on checmical structures.

The target proteins have to be represented by a binary vector to indicate whether the given protein is annotated as the drug target.

Both chemical properties and target profile must represent the same set of drugs, and prepare in a R list data, such as drug_tst in example.data.

3. Running DeepDR

The main function of DeepDR is DeepDRpredictor.R. Get your input files prepared, and run it like this:

Usage example:

predictions<-DeepDRpredictor(cell_tst,drug_tst) 

The output is the porbability of




