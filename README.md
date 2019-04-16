# DeepDRv1
DeepDRv1 includes the preditor (DeepDRpredictor.R) for prediction any given cell-drug relationship, DeepDRc.R for generating the best drug repurposing combination for cancer patient treatment.
It also includes procedure for generating the patients clone information from allele fraction data that came from FIREHOSE (FIREHOSE-tumor-vaf-processing).
Before using DeepDRv1, the example test data cell_tst, drug_tst, and mut_data.RData have to be obtained at /testdata.

#To get the prediction score (DeepDR score) for tesing cell-drug pairs, run the following in R

DeepDRpredictor<-DeepDRpredictor(cell_tst,drug_tst,"mutation")

#To get the best drug repurposing combination following the R code of DeepDRc.R by using colon-scMut-prediction.RData.

final.DC<-final_rDC(patient1_prediction_compare)

#To get the patients clone information, you have to install sciclone R package first

sciclone is available at https://github.com/genome/sciclone
