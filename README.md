# DeepDRK
**Deep** learning of **D**rug **R**esponse using **k**ernel-based data integration.

# Status

Active development

# Introduction

Early prediction of therapeutic response of cancer patients is a critical step in precision oncology. Among various obstacles hindering clinical translation, lacking effective methods for multimodal and multi-source data integration has become a bottleneck. DeepDRK provides a systematic way to predict drug response of personalized cancer cells via kernel-based data integration of pharmacogenomics, transcriptomics, epigenomics, chemical properties, and previously reported protein-compound interaction collected from different resources.

# Usage

1. Installation

Prerequisites: (1) R is properly installed; (2) Rscript is available in your system path ($PATH); (3) git (2.21.1)

Installation: git clone https://github.com/wangyc82/DeepDRv1

If the data (combination-data.RData) is incomplete, download it by clicking this file.

Dependencies: 

R;

Readr1.3.1 and all its dependencies;

h2o package (h2o_3.32.0.1.tgz) and its dependencies;

Oracle JDK (Java SE Development Kit 11.0.9).

Testing of successful installation by running the following in R

library(readr)

library(h2o)

h2o.init()

and get something like the following

![DeepDRKgithub1](https://user-images.githubusercontent.com/36029032/103272576-341d8900-49f8-11eb-9164-99cf3174aa95.png)



2. Preparation of the input files

Six csv files are needed to represent the multi-omic profile of cancer cells, i.e., single nucleotide variant and small INDELs (mutation.csv), copy number alteration (CN.csv), DNA methylation (Methy.csv), gene expression (Exp.csv), compound chemical properties (chem.csv), and known drug targets (DT.csv). 


The input file representing the genomic mutations of the cancer cells is a data matrix. The rows of the matrix are cancer cells, the columns are genes, and the elements of the matrix are variables of 0 or 1. 1 means the target gene is mutated in the cell, and 0 means no mutation in the cell. The number of column is the number of genes, and they are separated by tab delimiter.
 
mutation.csv

<img width="432" alt="Screenshot 2563-12-29 at 5 09 13 PM" src="https://user-images.githubusercontent.com/36029032/103272932-0c7af080-49f9-11eb-85c8-6c922dd6cefd.png">

The input file representing the copy number alterations, the status of DNA methylation, and the gene expression of the cancer cells are data matrices (CN.csv, methylation.csv, expression.csv). The rows of these matrices are cancer cells, the columns are genes, and the elements of the matrices are variables of float, integer, and float, respectively, representing the copy number, methylation, and expression of gene in the cell. The number of column of these matrices are the number of genes, and columns are separated by tab delimiter.

expression.csv

<img width="444" alt="Screenshot 2563-12-29 at 5 14 38 PM" src="https://user-images.githubusercontent.com/36029032/103273070-5ebc1180-49f9-11eb-9ac7-9174c71529d9.png">

The input file representing the chemical properties of the cancer drugs is a data matrix (chem.csv). The rows of the matrix are cancer drugs, the columns are chemical properties, and the elements of the matrix are variables of float. The number of column is the number of chemical properties, and the columns are separated by tab delimiter. To get chemical properties of cancer drugs, the chemical structure sdf file has to be prepared in advance, and upload it into “StarVue” (StarVue-macinstall-1.4.dmg) software to get the chemical molecular descriptors (QuaSAR-Descriptor in the Molecular Operating Environment (MOE)).

chem.csv

<img width="965" alt="Screenshot 2563-12-29 at 5 09 42 PM" src="https://user-images.githubusercontent.com/36029032/103273098-74c9d200-49f9-11eb-819a-fe6d82ceea6f.png">

The input file representing the known target proteins of the cancer drugs is a data matrix (DT.csv). The rows of the matrix are cancer drugs, the columns are target proteins, and the elements of the matrix are variables of 0 or 1. 1 means the target protein is a binding protein for a given drug, and 0 means not binding. The number of column is the number of proteins, and the columns are separated by tab delimiter.

DT.csv

<img width="432" alt="Screenshot 2563-12-29 at 5 09 13 PM" src="https://user-images.githubusercontent.com/36029032/103273223-b6f31380-49f9-11eb-8b91-70fa23e340fd.png">

The example input files can be found in Github repository data folder.

3. Running DeepDR

The main function of DeepDR is DeepDRpredictor.R. Get your input files prepared, and run it like this:

Usage example:

# loading the input files:

cell_tst<-list()

library(readr)

mutation <- read_csv("~/DeepDRv1/data/mutation.csv");A<-data.matrix(mutation[,-1]);rownames(A)<-mutation$X1;cell_tst[[1]]<-A

CN <- read_csv("~/DeepDRv1/data/CN.csv");A<-data.matrix(CN[,-1]);rownames(A)<- CN$X1;cell_tst[[2]]<-A

Methy <- read_csv("~/DeepDRv1/data/Methy.csv");A<-data.matrix(Methy[,-1]);rownames(A)<- Methy$X1;cell_tst[[3]]<-A

Exp <- read_csv("~/DeepDRv1/data/Exp.csv");A<-data.matrix(Exp[,-1]);rownames(A)<- Exp$X1;cell_tst[[4]]<-A

drug_tst<-list()

chem <- read_csv("~/DeepDRv1/data/chem.csv");A<-data.matrix(chem[,-1]);rownames(A)<- chem$X1;drug_tst[[1]]<-A

DT <- read_csv("~/DeepDRv1/data/DT.csv");A<-data.matrix(DT[,-1]);rownames(A)<- DT$X1;drug_tst[[2]]<-A

source('~/DeepDRv1/DeepDRpredictor.R')

predictions<-DeepDRpredictor(cell_tst,drug_tst)

<img width="348" alt="Screenshot 2563-12-29 at 11 49 02 AM" src="https://user-images.githubusercontent.com/36029032/103273130-8c08bf80-49f9-11eb-834a-d53cded05b17.png">

The true relationships between test cells and drugs are

<img width="324" alt="Screenshot 2563-12-29 at 11 52 23 AM" src="https://user-images.githubusercontent.com/36029032/103273169-9fb42600-49f9-11eb-98d5-1fbccb7a7097.png">


The predictors for the extender DeepDRK model dealing with the missing features (DeepDRpredictor.e) was also included. Here is the example showing how to use it:

suppose the mutation, methylation and target proteins are missing

cell_tst<-list()

library(readr)

CN <- read_csv("~/DeepDRv1/data/CN.csv");A<-data.matrix(CN[,-1]);rownames(A)<- CN$X1;cell_tst[[2]]<-A

Exp <- read_csv("~/DeepDRv1/data/Exp.csv");A<-data.matrix(Exp[,-1]);rownames(A)<- Exp$X1;cell_tst[[4]]<-A

drug_tst<-list()

chem <- read_csv("~/DeepDRv1/data/chem.csv");A<-data.matrix(chem[,-1]);rownames(A)<- chem$X1;drug_tst[[1]]<-A

drug_tst[[2]]<-matrix()

missCtype=c(1,3)

missDtype=2

source('~/DeepDRv1/DeepDRpredictor.e.R')

predictions<-DeepDRpredictor.e(cell_tst,drug_tst,missCtype,missDtype)

<img width="346" alt="Screenshot 2563-12-29 at 4 51 48 PM" src="https://user-images.githubusercontent.com/36029032/103273243-cbcfa700-49f9-11eb-9b86-79c91c8c6ff1.png">

# Contact

For technical issues please send an email to ycwang@nwipb.cas.cn.
