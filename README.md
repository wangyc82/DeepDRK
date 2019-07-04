# DeepDRv1
DeepDRv1 includes the preditor (DeepDRpredictor.R) for prediction any given cell-drug relationship, DeepDRc.R for generating the best drug repurposing combination for cancer patient treatment.

It also includes procedure for generating the patients clone information from allele fraction data that came from FIREHOSE (FIREHOSE-tumor-vaf-processing).

Before using DeepDRv1, the example test data cell_tst, drug_tst, and mut_data.RData have to be obtained at /testdata.

To get the prediction score (DeepDR score) for tesing cell-drug pairs with the monoclonal model, run the following in R

DeepDRpredictor<-DeepDRpredictor(cell_tst,drug_tst,"mutation")

To get the toxicity for test drug, using toxicity.R by the follwoing:
tox<-toxicity(drug_known,drug_test,AUCDR)

#AUCDR is the drug response matrix in training set with column as drug

#drug_known is the drug profile in training set

#drug_test is the drug profile in test set

To get the best drug repurposing combination with the poloclonal model running the following the R code with DeepDRc.R

final.DC<-DeepDRc(mono_score,cutoff,drug_tox)

#the mono_score is the monoclonal model prediction score with row as the drug

#drug_tox is the toxxicity for test drug

#cutoff is for distinguish sensitives from resistants

To get the patients clone information, you have to install sciclone R package first

sciclone is available at https://github.com/genome/sciclone
