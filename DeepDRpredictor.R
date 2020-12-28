# this procedure is for prediction of the relationship between cancer cells/patients and anti-cancer drugs
# it could be applied on mutation, copy number variation, DNA methylation, gene expression profiles for cell/patients, and drug chemical properties and target proteins
# drug chemical properties are generated through QuaSAR-Descriptor in the Molecular Operating Environment (MOE) based on drug's chemical structure
# the cell/patient data is deposited as a list in cell_tst containg mutation, copy number variation, DNA methylation, and gene expression profiles for common cells/patients
# the drug data is deposited as a list in drug_tst containg chemical properties, and target proteins for common drugs
DeepDRpredictor<-function(cell_tst,drug_tst)
{
  load("~/DeepDRv1/combination-data.RData")
  
  #preparing the inputs for input for deep learning method
  sim_cell<-cbind(sim_mut,sim_CN,sim_methy,sim_exp)
  sim_drug<-cbind(sim_comp,sim_DT)
  
  A<-which(AUC_matS_comb==1,arr.ind = TRUE)
  B<-which(AUC_matS_comb==-1,arr.ind = TRUE)
  Xprn<-cbind(sim_cell[A[,1],],sim_drug[A[,2],])
  Xnrn<-cbind(sim_cell[B[,1],],sim_drug[B[,2],])
  Xrn<-rbind(Xprn,Xnrn)
  Yrn<-c(rep(0,nrow(A)),rep(1,nrow(B)))
  data_trn<-data.frame(Xrn,Yrn)
  library(h2o)
  h2o.init()
  cat("convert training data \n")
  train.hex <- as.h2o(data_trn)
  #prepare for testing
  #make sure they have same cells and drugs
  
  cell_trs_mut<-cell_tst[[1]]
  cell_trs_CN<-cell_tst[[2]]
  cell_trs_methy<-cell_tst[[3]]
  cell_trs_exp<-cell_tst[[4]]
  
  Csim_mut<-exp(-as.matrix(dist(rbind(cell_mut[,intersect(colnames(cell_mut),colnames(cell_trs_mut))],cell_trs_mut[,intersect(colnames(cell_mut),colnames(cell_trs_mut))]),method = "binary")))
  Csim_CN<-exp(-0.001*as.matrix(dist(rbind(cell_CN[,intersect(colnames(cell_CN),colnames(cell_trs_CN))],cell_trs_CN[,intersect(colnames(cell_CN),colnames(cell_trs_CN))]))))
  Csim_methy<-exp(-0.001*as.matrix(dist(rbind(cell_methy[,intersect(colnames(cell_methy),colnames(cell_trs_methy))],cell_trs_methy[,intersect(colnames(cell_methy),colnames(cell_trs_methy))]))))
  Csim_exp<-exp(-0.001*as.matrix(dist(rbind(cell_exp[,match(intersect(colnames(cell_exp),colnames(cell_trs_exp)),colnames(cell_exp))],cell_trs_exp[,match(intersect(colnames(cell_exp),colnames(cell_trs_exp)),colnames(cell_trs_exp))]))))
  
  Csim_tst<-cbind(Csim_mut[rownames(cell_trs_mut),rownames(cell_mut)],Csim_CN[rownames(cell_trs_CN),rownames(cell_CN)],Csim_methy[rownames(cell_trs_methy),rownames(cell_methy)],Csim_exp[rownames(cell_trs_exp),rownames(cell_exp)])
  
  drug_tst_comp<-drug_tst[[1]]
  drug_tst_DT<-drug_tst[[2]]
  
  Dsim_comp<-exp(-0.001*as.matrix(dist(rbind(drug_comp[,intersect(colnames(drug_comp),colnames(drug_tst_comp))],drug_tst_comp[,intersect(colnames(drug_comp),colnames(drug_tst_comp))]))))
  Dsim_DT<-exp(-as.matrix(dist(rbind(drug_DT[,intersect(colnames(drug_DT),colnames(drug_tst_DT))],drug_tst_DT[,intersect(colnames(drug_DT),colnames(drug_tst_DT))]),method = "binary")))
  
  Dsim_tst<-cbind(Dsim_comp[rownames(drug_tst_comp),rownames(drug_comp)],Dsim_comp[rownames(drug_tst_DT),rownames(drug_DT)])
  
  Ast<-which(matrix(1,nrow(cell_trs_mut),nrow(drug_tst_comp))==1,arr.ind = TRUE)
  Xst<-cbind(Csim_tst[Ast[,1],],Dsim_tst[Ast[,2],])
  Yst<-sample(c(0,1), nrow(Xst), replace = TRUE)
  data_tst<-data.frame(Xst,Yst)
  colnames(data_tst)<-colnames(data_trn)
  cat("convert testing data \n")
  test.hex <- as.h2o(data_tst)
  cat("bulid the training model \n")
  model=h2o.deeplearning(x = 1:ncol(Xrn), y = ncol(Xrn)+1, training_frame = train.hex, validation=test.hex,hidden=c(200,200), epochs=10, activation="Tanh")
  cat("performing the prediction \n")
  model_prediction<-h2o.predict(model, test.hex)
  prob_test<-as.data.frame(model_prediction)[,1]
  return(prob_test)
  h2o.shutdown()
}
  
