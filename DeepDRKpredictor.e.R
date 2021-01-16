# This program is to handle task with missing features  
# In case the mutation, methylation and target proteins are missing, tt requires some following input to characterize a cancer cell line: 
# 1. Copy number variation (CN.csv);
# 2. Gene expression profile (Exp.csv);
# Meanwhile, it requires the following input to characterize a drug: 
# 3. Drug chemical properties, generated through QuaSAR-Descriptor in the Molecular Operating Environment (MOE) based on drug's chemical structure (chem.csv); 
 
# Example to run:
# cell_tst<-list()

# library(readr)

# CN <- read_csv("~/DeepDRK/data/CN.csv")
# A<-data.matrix(CN[,-1])
# rownames(A)<- CN$X1
# cell_tst[[2]]<-A

# Exp <- read_csv("~/DeepDRK/data/Exp.csv")
# A<-data.matrix(Exp[,-1])
# rownames(A)<- Exp$X1
# cell_tst[[4]]<-A

# drug_tst<-list()

# chem <- read_csv("~/DeepDRK/data/chem.csv")
# A<-data.matrix(chem[,-1])
# rownames(A)<- chem$X1
# drug_tst[[1]]<-A
# drug_tst[[2]]<-matrix()

# missCtype=c(1,3)
# missDtype=2
# source('~/DeepDRK/DeepDRKpredictor.e.R')
# predictions<-DeepDRKpredictor.e(cell_tst,drug_tst,missCtype,missDtype)
DeepDRKpredictor.e<-function(cell_tst,drug_tst,missCtype,missDtype)
{
  #load("~/DeepDRv1/combination-data.RData")
  
  #preparethe inputs for deep learning method
  sim_cell_ls<-list()
  sim_cell_ls[[1]]<-sim_mut
  sim_cell_ls[[2]]<-sim_CN
  sim_cell_ls[[3]]<-sim_methy
  sim_cell_ls[[4]]<-sim_exp
  sim_cell_ls1<-sim_cell_ls[-missCtype]
  
  sim_cell<-NULL;for (i in 1:length(sim_cell_ls1)) {sim_cell<-cbind(sim_cell,sim_cell_ls1[[i]])}
  
  sim_drug_ls<-list()
  sim_drug_ls[[1]]<-sim_comp
  sim_drug_ls[[2]]<-sim_DT
  sim_drug<-sim_drug_ls[[-missDtype]]
  
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
  test_cell<-unique(c(rownames(cell_trs_mut),rownames(cell_trs_CN),rownames(cell_trs_methy),rownames(cell_trs_exp)))
  
  Csim_mut<-exp(-as.matrix(dist(rbind(cell_mut[,intersect(colnames(cell_mut),colnames(cell_trs_mut))],cell_trs_mut[,intersect(colnames(cell_mut),colnames(cell_trs_mut))]),method = "binary")))
  Csim_CN<-exp(-0.001*as.matrix(dist(rbind(cell_CN[,intersect(colnames(cell_CN),colnames(cell_trs_CN))],cell_trs_CN[,intersect(colnames(cell_CN),colnames(cell_trs_CN))]))))
  Csim_methy<-exp(-0.1*as.matrix(dist(rbind(cell_methy[,intersect(colnames(cell_methy),colnames(cell_trs_methy))],cell_trs_methy[,intersect(colnames(cell_methy),colnames(cell_trs_methy))]))))
  Csim_exp<-exp(-0.01*as.matrix(dist(rbind(cell_exp[,match(intersect(colnames(cell_exp),colnames(cell_trs_exp)),colnames(cell_exp))],cell_trs_exp[,match(intersect(colnames(cell_exp),colnames(cell_trs_exp)),colnames(cell_trs_exp))]))))
  
  Csim_ls<-list()
  Csim_ls[[1]]<-Csim_mut
  Csim_ls[[2]]<-Csim_CN
  Csim_ls[[3]]<-Csim_methy
  Csim_ls[[4]]<-Csim_exp
  Csim_ls1<-Csim_ls[-missCtype]
  Csim_tst<-NULL
  for (i in 1:length(Csim_ls1)) {
    MM<-Csim_ls1[[i]]
    Csim_tst<-cbind(Csim_tst,MM[test_cell,rownames(cell_mut)])
    rm(MM)
  }
  
  drug_tst_comp<-drug_tst[[1]]
  drug_tst_DT<-drug_tst[[2]]
  test_drug<-unique(c(rownames(drug_tst_comp),rownames(drug_tst_DT)))
  
  Dsim_comp<-exp(-0.001*as.matrix(dist(rbind(drug_comp[,intersect(colnames(drug_comp),colnames(drug_tst_comp))],drug_tst_comp[,intersect(colnames(drug_comp),colnames(drug_tst_comp))]))))
  Dsim_DT<-exp(-as.matrix(dist(rbind(drug_DT[,intersect(colnames(drug_DT),colnames(drug_tst_DT))],drug_tst_DT[,intersect(colnames(drug_DT),colnames(drug_tst_DT))]),method = "binary")))
  
  Dsim_ls<-list()
  Dsim_ls[[1]]<-Dsim_comp[test_drug,rownames(drug_comp)]
  Dsim_ls[[2]]<-Dsim_DT[test_drug,rownames(drug_comp)]
  Dsim_tst<-Dsim_ls[[-missDtype]]
   
  Ast<-which(matrix(1,length(test_cell),length(test_drug))==1,arr.ind = TRUE)
  Xst<-cbind(Csim_tst[Ast[,1],],Dsim_tst[Ast[,2],])
  Yst<-sample(c(0,1), nrow(Xst), replace = TRUE)
  data_tst<-data.frame(Xst,Yst)
  colnames(data_tst)<-colnames(data_trn)
  cat("convert testing data \n")
  test.hex <- as.h2o(data_tst)
  cat("bulid the training model \n")
  model=h2o.deeplearning(x = 1:ncol(Xrn), y = ncol(Xrn)+1, training_frame = train.hex, validation=test.hex,hidden=c(200,200), epochs=10, activation="Tanh",seed=1, reproducible=T)
  cat("performing the prediction \n")
  model_prediction<-h2o.predict(model, test.hex)
  prob_test<-as.data.frame(model_prediction)[,1]
  S<-data.frame(cell=test_cell[Ast[,1]],drug=test_drug[Ast[,2]],prob=prob_test)
  return(S)
  h2o.shutdown()
}

