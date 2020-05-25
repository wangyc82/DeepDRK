# DeepDRv1
DeepDRv1 includes the preditor (DeepDRpredictor.R) for prediction any given cell-drug relationship, DeepDRc.R for generating the best drug repurposing combination for cancer patient treatment.

It also includes procedure for generating the patients clone information from allele fraction data that came from FIREHOSE (FIREHOSE-tumor-vaf-processing).

To get the prediction score (DeepDR score) for tesing cell-drug pairs with the monoclonal model, run the following in R

DeepDRpredictor<-DeepDRpredictor(cell_tst,drug_tst,"mutation")

#cell_tst is test cell mutation profile

#drug_tst is test drug chemical properties profile

To get the toxicity for test drug, using toxicity.R by the follwoing:
tox<-toxicity(drug_known,drug_test,AUCDR)

#AUCDR is the drug response matrix in training set with column as drug

#drug_known is the drug profile in training set

#drug_test is the drug profile in test set
