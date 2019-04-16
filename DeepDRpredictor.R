#this this procedure is for prediction of the relationship between cancer cells/patients and anti-cancer drugs
#this is based on CTRP and GDSC non-redundent data
#it could be applied either for expression or mutation cell/patients profile, and drug chemical property
#drug chemical property is generated through QuaSAR-Descriptor in the Molecular Operating Environment (MOE)
#categary is indicater for either expression or mutation
# make sure you cell_profile and drug profile have the same feature (same gene set for cell_profile, same descriptor for drug_profile) with the training data
# the cell/patient profile is as cell_tst
# the drug profile is as drug_tst
DeepDRpredictor<-function(cell_tst,drug_tst,categary)
{
  if (categary=="mutation") {load("~/Documents/DeepDR/version at HongKong/DeepDR-sofware/mut-data.RData")} else {load("~/Documents/DeepDR/version at HongKong/DeepDR-sofware/exp-data.RData")}
  #preparing for training
  A<-which(AUCmat_comb_S==1,arr.ind = TRUE)
  B<-which(AUCmat_comb_S==-1,arr.ind = TRUE)
  Xprn<-cbind(sim_cell[A[,1],],sim_comp[A[,2],])
  Xnrn<-cbind(sim_cell[B[,1],],sim_comp[B[,2],])
  Xrn<-rbind(Xprn,Xnrn)
  Yrn<-c(rep("a",nrow(A)),rep("b",nrow(B)))
  data_trn<-data.frame(Xrn,Yrn)
  library(h2o)
  h2o.init()
  cat("convert training data \n")
  train.hex <- as.h2o(data_trn)
  #prepare for testing
  if (categary=="mutation") {Csim<-exp(-as.matrix(dist(rbind(cell_pro,cell_tst),method = "binary")))} else {Csim<-exp(-0.001*as.matrix(dist(rbind(cell_pro,cell_tst))))}
  Csim_tst<-Csim[rownames(cell_tst),rownames(cell_pro)]
  Dsim<-exp(-0.001*as.matrix(dist(rbind(drug_comp,drug_tst))))
  Dsim_tst<-Dsim[rownames(drug_tst),rownames(drug_comp)]
  Ast<-which(matrix(1,nrow(cell_tst),nrow(drug_tst))==1,arr.ind = TRUE)
  Xst<-cbind(Csim_tst[Ast[,1],],Dsim_tst[Ast[,2],])
  Yst<-sample(c("a","b"), nrow(Xst), replace = TRUE)
  data_tst<-data.frame(Xst,Yst)
  colnames(data_tst)<-colnames(data_trn)
  cat("convert testing data \n")
  test.hex <- as.h2o(data_tst)
  cat("bulid the training model \n")
  model=h2o.deeplearning(x = 1:ncol(Xrn), y = ncol(Xrn)+1, training_frame = train.hex, validation=test.hex,hidden=c(200,200), epochs=10, activation="Tanh")
  cat("performing the prediction \n")
  model_prediction<-h2o.predict(model, test.hex)
  prob_test<-as.data.frame(model_prediction)[,3]
  return(prob_test)
  h2o.shutdown()
}
