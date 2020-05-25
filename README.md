# DeepDRv1
DeepDRv1 includes the preditor (DeepDRpredictor.R) for prediction any given cell-drug relationship.

To get the prediction score (DeepDR score) for tesing cell-drug pairs with the monoclonal model, run the following in R

Usage example:

DeepDRpredictor<-DeepDRpredictor(cell_tst,drug_tst,categary) 

#cell_tst is test cell profile, could be one of mutation, copy number, methylation,expression, or integration

#drug_tst is test drug chemical properties profile

#categary is the data type for cancer genomics, could be one of mutation, copy number, methylation,expression, or integration

This procedure is used CTRP and GDSC non-redundant data (including 10,754 sensitive and 10,607 resistant cell-drug pairs) as the training data

#drug_known is the drug profile in training set

#cell_known is the cell profile in training set

#AUCmat_comb_S is the digitalized drug sensitivities 
