# DeepDRv1
DeepDRv1 includes the preditor (DeepDRpredictor.R) for prediction any given cell-drug relationship.

To get the prediction score (DeepDR score) for tesing cell-drug pairs, just run the DeepDRpredictor.R.

First you have to prepare the omics profile for cancer cells, including mutation at gene level, copy number at gene level, DNA methylation at gene level, and gene expression.

All these omics profiles must represent same set of cancer cells, that is, you have to prepare these omics data for same set of cancer cells, and let them is a R list data, such as cell_tst.

Then you have to prepare the profiles for drugs, including chemical properties generated from chemical structures, and target proteins.

The target proteins have to be represented by a binary vector to indicate whether the given protein is annotated as the target of this drug.

All chemical properties and target profile must represent the same set of drugs, and prepare in a R list data, such as drug_tst.

To use the predictor (DeepDRpredictor.R), the benchmark data for training have to be downloaded in advance. 

The benchmark data for training is combination-data.RData.

combination-data.RData includes the cancer-drug associations generated from GDSC and CTRP non-redundent data. 

combination-data.RData includes cancer omics profiles and drug chemical property and target protein profiles.

combination-data.RData includes the cancer and drug similarity matrices based on cancer omics and drug chemical properties and target proteins.

Usage example:

DeepDRpredictor<-DeepDRpredictor(cell_tst,drug_tst) 




