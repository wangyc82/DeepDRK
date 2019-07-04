#this is procedure to calculate the toxicity of drugs
# AUCDR is the drug response matrix in training set with column as drug
# drug_known is the drug profile in training set
# drug_test is the drug profile in test set
toxicity<-function(drug_known,drug_test,AUCDR)
{
  DD<-rbind(drug_known[,intersect(colnames(drug_known),colnames(drug_test))],drug_test[,intersect(colnames(drug_known),colnames(drug_test))])
  # construction the toxicity for known drug in training set
  Ti<-colSums(AUCDR)/nrow(AUCDR)
  # construction the toxicity for unknown drug in CTRP
  lab<-match(rownames(drug_test),rownames(drug_known))
  Ti.tst<-as.numeric()
  ll1<-which(lab!=0)
  ll2<-seq(1,length(lab),1)[-ll1]
  for (i in 1:length(ll1)) {Ti.tst[ll1[i]]<-Ti[lab[ll1[i]]]}
  sim_drug<-exp(-0.001*as.matrix(dist(DD)))
  sim_drug_test<-sim_drug[(nrow(drug_known)+1):nrow(DD),1:nrow(drug_known)]
  for (i in 1:length(ll2)) {Ti.tst[ll2[i]]<-sim_drug_test[ll2[i],]%*%Ti/length(Ti)}
  names(Ti.tst)<-rownames(drug_test)
  toxicity<-Ti.tst
}