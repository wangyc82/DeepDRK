# DeepDR
**Deep** learning **D**rug **R**esponses.

# Status

Active development

# Introduction

Accurately predicting the response of a cancer patient to a therapeutic agent is core goal of precision oncology. Among various obstacles hindering clinical translation, the lack of effective multimodal and multi-source data integration methods has become a bottleneck. DeepDR provides a systematic way to predict drug response from different resources, diverse cancer types, and various chemical compounds, based on kernel-based integration of multi-source data from genomics,
transcriptomics, epigenomics, and chemical properties and previously reported target proteins of compounds.

# Usage

1. Installation

Before installing DeepDR please make sure you have installed R, and Rscript is available in your system path ($PATH).

git clone https://github.com/wangyc82/DeepDRv1

2. Preparing the input files

The input files includes the omics profiles for cancer cells, including mutation at gene level, copy number at gene level, DNA methylation at gene level, and gene expression, and the drug properties, including chemical properties generated from chemical structures, and previously reported target proteins.

All these omics profiles must represent the same set of cancer cells, and they are prepared into R list format, such as cell_tst in example.data.

The chemical properties are the molecular properties based on checmical structures, and the target proteins have to be represented by a binary vector to indicate whether the given protein is annotated as the drug target.

Both chemical properties and target profile must represent the same set of drugs, and have to be prepared into R list format, such as drug_tst in example.data.

3. Running DeepDR

The main function of DeepDR is DeepDRpredictor.R. Get your input files prepared, and run it like this:

Usage example:

predictions<-DeepDRpredictor(cell_tst,drug_tst) 

The output describes the probability that the test cell is sensitive to the test drug.

# Dependencies

The model is trained by using h2o package in R. The dependency requirements are automatically solved while running the program.

# Contact

For technical issues please send an email to ycwang@nwipb.cas.cn or jgwang@ust.hk.






