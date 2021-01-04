<img src="https://user-images.githubusercontent.com/36029032/103525552-60cf1600-4eba-11eb-910e-b4a0ffb3e018.png" width="200">

**Deep** learning of **D**rug **R**esponse using **k**ernel-based data integration.

# Status

Active development

# Introduction

Early prediction of therapeutic response of cancer patients is a critical step in precision oncology. Among various obstacles hindering clinical translation, lacking effective methods for multimodal and multi-source data integration has become a bottleneck. DeepDRK provides a systematic way to predict drug response of personalized cancer cells via kernel-based data integration of pharmacogenomics, transcriptomics, epigenomics, chemical properties, and previously reported protein-compound interaction collected from different resources.

# Usage

1. Installation

   Prerequisites of DeepDRK including: 

   - R is properly installed; 

   - Rscript is available in your system path ($PATH);

   - git (2.21.1)

    Installation: 

    - step 1: git clone https://github.com/wangyc82/DeepDRK;

    - step 2: download the combined training data (combination_data.RData) from https://wanglab.shinyapps.io/DeepDRK/, and put it in the DeepDRK folder.

    Dependencies of DeepDRK including 

    - Readr1.3.1 and all its dependencies;

    - Oracle JDK (Java SE Development Kit 11.0.9);

    - h2o package (h2o_3.32.0.1.tgz) and its dependencies.

    Testing of successful installation by running the following commands in R:
     > library(readr)
     > h2o.init()


2. Preparation of the input files

In total, six csv files are needed to represent the multi-omic profile of cancer cells, i.e., single nucleotide variant and small INDELs (mutation.csv), copy number alteration (CN.csv), DNA methylation (methylation.csv), gene expression (expression.csv), compound chemical properties (chem.csv), and known drug targets (DT.csv). 

The input file mutation.csv that includes genomic mutations in cancer cells is a data matrix, with each row representing a cancer cell line and each column representing the genotype of a gene. The value of this matrix is binary, with 1 indicating mutated while 0 for wild type.
 
mutation.csv

              A1BG A1CF  A2M …
    201T       0     0    0
    22RV1      0     1    0
    42-MG-BA   0     0    1
      .
      .
      .

Similarly, the input files that contain the copy number alteration data (CN.csv), the status of DNA methylation (methylation.csv), and the gene expression of the cancer cells (expression.csv) are also data matrices with a row representing a cancer cell line and a column representing a gene. The elements of these matrices are respectively integers for gene copy numbers, and float numbers for level of gene methylation and expression.

expression.csv

               A1BG   A1CF    A2M …
    201T      3.162   2.919   3.379
    22RV1     3.531   6.336   5.331
    42-MG-BA  6.002   3.137   3.237
      .
      .
      .

The input file chem.csv that describes the chemical properties of the cancer drugs is a matrix with each row representing one cancer drug and each column representing one feature to describe drug’s chemical properties. The descriptors of the chemical properties of a cancer drug were inferred from its chemical structure. Particularly, to describe a drug, such as Erlotinib, we will need to first download the sdf file of this drug from PubChem, and then upload the chemical structure into “StarVue” (StarVue-macinstall-1.4.dmg) software to extract the 2D Molecular Operating Environment (MOE)) descriptors, including physical properties, atom counts, and bond counts.

chem.csv

              PUBCHEM_MOLECULAR_WEIGHT   PUBCHEM_EXACT_MASS    PUBCHEM_CACTVS_TPSA …
    Erlotinib            393.4                   393.2                  74.4
    Rapamycin            917.2                   913.6                  195.0
    Sunitinib            398.5                   398.2                  77.2
       .
       .
       .


The input file DT.csv includes the known targeting proteins of the cancer drugs. Each row represents one cancer drug, and each column represents a target protein. “1” indicates a potential drug-gene interaction reported in DrugBank or KEGG.

DT.csv

                         EGFR                     KIT                   PDGRA …
    Erlotinib            393.4                   393.2                  74.4
    Rapamycin            917.2                   913.6                  195.0
    Sunitinib            398.5                   398.2                  77.2
       .
       .
       .

The example of all input files can be found in the “data” folder of the Github repository.

3. Running DeepDRK

The main function of DeepDRK is DeepDRKpredictor.R. Get your input files prepared, and run it like this:

Usage example:

    > cell_tst<-list()
    > library(readr)
    > mutation <- read_csv("~/DeepDRK/data/mutation.csv");A<-data.matrix(mutation[,-1]);rownames(A)<-mutation$X1;cell_tst[[1]]<-A
    > CN <- read_csv("~/DeepDRK/data/CN.csv");A<-data.matrix(CN[,-1]);rownames(A)<- CN$X1;cell_tst[[2]]<-A
    > Methy <- read_csv("~/DeepDRK/data/methylation.csv");A<-data.matrix(Methy[,-1]);rownames(A)<- Methy$X1;cell_tst[[3]]<-A
    > Exp <- read_csv("~/DeepDRK/data/expression.csv");A<-data.matrix(Exp[,-1]);rownames(A)<- Exp$X1;cell_tst[[4]]<-A
    > drug_tst<-list()
    > chem <- read_csv("~/DeepDRK/data/chem.csv");A<-data.matrix(chem[,-1]);rownames(A)<- chem$X1;drug_tst[[1]]<-A
    > DT <- read_csv("~/DeepDRK/data/DT.csv");A<-data.matrix(DT[,-1]);rownames(A)<- DT$X1;drug_tst[[2]]<-A
    > load("~/DeepDRK/combination_data.RData") #load the training RData
    > source('~/DeepDRK/DeepDRKpredictor.R')
    > predictions<-DeepDRKpredictor(cell_tst,drug_tst)
          cell       drug       prob
     1    697      Imatinib    0.881
     2  A3-KAW     Imatinib    0.819
     3    697   Gemcitabline   0.122
     4  A3-KAW  Gemcitabline   0.158
     > h2o.shutdown() # shut down the h2o
     Are you sure you want to shutdown the H2O instance running at http://localhost:54321/ (Y/N)? y
     TRUE
     
As shown in the following figure of the experimental data, we observed that cell lines 697 and A3-KAW tend to be sensitive to Imatinib and these two cell lines were prune to be resistant to Gemcitabine, consistent with the prediction from DeepDRK. 

![example-test-AUCDR](https://user-images.githubusercontent.com/36029032/103406015-1dac3480-4b94-11eb-8981-31293cd1d231.png)

Moreover, DeepDRK could also handle task with missing features using the DeepDRKpredictor.e R function. Here is the example showing how to use it:

In case the mutation, methylation and target proteins are missing

    > cell_tst<-list()
    > library(readr)
    > CN <- read_csv("~/DeepDRK/data/CN.csv");A<-data.matrix(CN[,-1]);rownames(A)<- CN$X1;cell_tst[[2]]<-A
    > Exp <- read_csv("~/DeepDRK/data/expression.csv");A<-data.matrix(Exp[,-1]);rownames(A)<- Exp$X1;cell_tst[[4]]<-A
    > drug_tst<-list()
    > chem <- read_csv("~/DeepDRK/data/chem.csv");A<-data.matrix(chem[,-1]);rownames(A)<- chem$X1;drug_tst[[1]]<-A
    > drug_tst[[2]]<-matrix()
    > missCtype=c(1,3)
    > missDtype=2
    > load("~/DeepDRK/combination_data.RData") #load the training RData
    > source('~/DeepDRK/DeepDRKpredictor.e.R')
    > predictions<-DeepDRKpredictor.e(cell_tst,drug_tst,missCtype,missDtype)
          cell       drug       prob
     1    697      Imatinib    0.706
     2  A3-KAW     Imatinib    0.721
     3    697   Gemcitabline   0.127
     4  A3-KAW  Gemcitabline   0.125

    > h2o.shutdown() # shut down the h2o
    Are you sure you want to shutdown the H2O instance running at http://localhost:54321/ (Y/N)? y
    TRUE


# Contact

For technical issues please send an email to ycwang@nwipb.cas.cn.
