# This program is to predict potential drug response of a cancer cell line (cancer patient)
# It requires the following input to characterize a cancer cell line: 
# 1. Gene mutation (mutation.csv);
# 2. Copy number variation (CN.csv);
# 3.  DNA methylation (Methy.csv);
# 4. Gene expression profile (Exp.csv);
# Meanwhile, it requires the following input to characterize a drug: 
# 5. Drug chemical properties, generated through QuaSAR-Descriptor in the Molecular Operating Environment (MOE) based on drug's chemical structure (chem.csv); 
# 6. Known drug targets (DT.csv).
# 
# Example to run:
# cell_tst<-list()
# library(readr)
# mutation <- read_csv("~/DeepDRK/data/mutation.csv")
# A<-data.matrix(mutation[,-1])
# rownames(A)<-mutation$X1
# cell_tst[[1]]<-A

# CN <- read_csv("~/DeepDRK/data/CN.csv")
# A<-data.matrix(CN[,-1])
# rownames(A)<- CN$X1
# cell_tst[[2]]<-A

# Methy <- read_csv("~/DeepDRK/data/Methy.csv")
# A<-data.matrix(Methy[,-1])
# rownames(A)<- Methy$X1
# cell_tst[[3]]<-A

# Exp <- read_csv("~/DeepDRK/data/Exp.csv")
# A<-data.matrix(Exp[,-1])
# rownames(A)<- Exp$X1
# cell_tst[[4]]<-A

# drug_tst<-list()
# chem <- read_csv("~/DeepDRK/data/chem.csv")
# A<-data.matrix(chem[,-1])
# rownames(A)<- chem$X1
# drug_tst[[1]]<-A

# DT <- read_csv("~/DeepDRK/data/DT.csv")
# A<-data.matrix(DT[,-1])
# rownames(A)<- DT$X1
# drug_tst[[2]]<-A

# source('~/DeepDRv1/DeepDRKpredictor.R')
# predictions<-DeepDRKpredictor(cell_tst,drug_tst)
DeepDRKpredictor<-function(cell_tst,drug_tst)
{
  #load("~/DeepDRK/combination-data.RData")
  
  # prepare the inputs for deep learning method
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
  
  #prepare the test data
  
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
  S<-data.frame(cell=rownames(cell_trs_exp)[Ast[,1]],drug=rownames(drug_tst_DT)[Ast[,2]],prob=prob_test)
  return(S)
  h2o.shutdown()
}
  
