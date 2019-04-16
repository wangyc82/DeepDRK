#this procedure is to get tumor vaf data from FIREHOSE
#this vaf data is for geting combination of clones for each patient based on sciClone
# only get the vaf data from BRCA, ESCA,SARC, TGCT,and UCEC
vaf_column<-c("Chromosome","Start_position","i_tumor_ref_reads","i_tumors_var_reads","i_tumor_vaf")
TCGA_barcode<-NULL
library(sciClone)
sc<-list()
MAXCLUSTER=10
plotIntermediateResults = 0
path<-"~/Documents/DeepDRpaper-version2/versionHongKong2018/sciClone/BRCA.Mutation_Packager_Oncotated_Calls.Level_3"
files <- dir(path, recursive=TRUE, full.names=TRUE, pattern="\\.txt$")

vaf_dat<-list()
for (i in 1:length(files)) {
  df<-read_delim(files[i],  "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)
  TCGA_barcode[i]<-unique(df$Tumor_Sample_Barcode)
  annot = df[,c(1,5,6)]
  df1<-df[,vaf_column]
  df2<-data.frame(df1)
  df3<-df2
  A<-unlist(gregexpr("\\|",df2$i_tumor_ref_reads))
  B<-unlist(gregexpr("\\|",df2$i_tumors_var_reads))
  C<-unlist(gregexpr("\\|",df2$i_tumor_vaf))
  if (length(which(A!=-1))==0) {
    df3$i_tumor_ref_reads<-df2$i_tumor_ref_reads
  } else {
    a<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_tumor_ref_reads,"\\|")[[x]][1]))
    df3$i_tumor_ref_reads<-as.numeric(a)
  }
  if (length(which(B!=-1))==0) {
    df3$i_tumors_var_reads<-df2$i_tumors_var_reads
  } else {
    b<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_tumors_var_reads,"\\|")[[x]][1]))
    df3$i_tumors_var_reads<-as.numeric(b)
  }
  if (length(which(C!=-1))==0) {
    df3$i_tumor_vaf<-df2$i_tumor_vaf
  } else {
    c<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_tumor_vaf,"\\|")[[x]][1]))
    df3$i_tumor_vaf<-as.numeric(c)
  }
  vaf_dat[[i]]<-data.frame(gene=df$Hugo_Symbol,df3)
  sc[[i]] <- sciClone(vafs=df3, sampleNames=unique(df$Tumor_Sample_Barcode), useSexChrs=FALSE, copyNumberCalls=NULL, copyNumberMargins=0.4, minimumDepth=0, verbose=1, doClusteringAlongMargins=F, plotIntermediateResults = plotIntermediateResults,maximumClusters = MAXCLUSTER,annotation = annot, clusterMethod = 'gaussian.bmm')
  rm(A,B,C)
  rm(df,df1,df2,df3,annot)
  cat(i,"\n")
  }
names(sc)<-TCGA_barcode
TCGA_patient<-substr(TCGA_barcode,1,12)
u_TCGApatient<-unique(TCGA_patient)
#get patient vaf data
vaf_dat_p<-lapply(1:length(u_TCGApatient),function(x) vaf_dat[which(TCGA_patient %in% u_TCGApatient[x])])
names(vaf_dat_p)<-u_TCGApatient
vaf_datp_gene<-lapply(1:length(vaf_dat_p),function(x) as.character(vaf_dat_p[[x]][[1]][,1]))
all_gene<-unique(unlist(vaf_datp_gene))
pat_vaf<-matrix(0,length(u_TCGApatient),length(all_gene))
rownames(pat_vaf)<-u_TCGApatient
colnames(pat_vaf)<-all_gene
for (i in 1:length(u_TCGApatient)) {pat_vaf[i,vaf_datp_gene[[i]]]<-as.numeric(vaf_dat_p[[i]][[1]][,6])}

# merge clone infomation for each patient
sc_p<-lapply(1:length(u_TCGApatient),function(x) sc[which(TCGA_patient %in% u_TCGApatient[x])])
names(sc_p)<-u_TCGApatient
vaf_dat_p<-lapply(1:length(u_TCGApatient),function(x) vaf_dat[which(TCGA_patient %in% u_TCGApatient[x])])
names(vaf_dat_p)<-u_TCGApatient
m<-NULL;for (i in 1:length(sc_p)) {m[i]<-length(sc_p[[i]][[1]])}
ig_id<-which(m==0) # can not get clone combinations because leass of ten samples in vaf dara
sc_p1<-sc_p[-ig_id]
vaf_dat_p1<-vaf_dat_p[-ig_id]
#get clone-gene mutation matrix for each patient
cp_mut_mat<-list()
vaf_datp1_gene<-list()
for (i in 1:length(vaf_dat_p1)) {
  A<-vaf_dat_p1[[i]][[1]]
  p<-which(A$Chromosome=="X" | A$Chromosome=="Y")
  if (length(p)==0) {B<-A} else {B<-A[-p,]}
  vaf_datp1_gene[[i]]<-B$gene
  C<-sc_p1[[i]][[1]]@clust$cluster.assignments
  M<-matrix(0,nrow(B),length(levels(as.factor(C))))
  rownames(M)<-B$gene
  for (k in 1:nrow(B)) {M[k,C[k]]<-1}
  cp_mut_mat[[i]]<-M
  rm(A,B,C,M)
  cat(i,"\n")
}
names(cp_mut_mat)<-names(sc_p1)

path<-"~/Documents/DeepDRpaper-version2/versionHongKong2018/sciClone/ESCA.Mutation_Packager_Oncotated_Calls.Level_3"
vaf_column<-c("Chromosome","Start_position","i_tumor_ref_reads","i_tumor_var_reads","i_tumor_VAF")
library(sciClone)
MAXCLUSTER=10
plotIntermediateResults = 0

TCGA_barcode<-NULL
sc<-list()
vaf_dat<-list()
for (i in 1:length(files)) {
  df<-read_delim(files[i],  "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)
  TCGA_barcode[i]<-unique(df$Tumor_Sample_Barcode)
  annot = df[,c(1,5,6)]
  df1<-df[,vaf_column]
  df2<-data.frame(df1)
  df3<-df2
  A<-unlist(gregexpr("\\|",df2$i_tumor_ref_reads))
  B<-unlist(gregexpr("\\|",df2$i_tumor_var_reads))
  C<-unlist(gregexpr("\\|",df2$i_tumor_VAF))
  if (length(which(A!=-1))==0) {
    df3$i_tumor_ref_reads<-df2$i_tumor_ref_reads
  } else {
    a<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_tumor_ref_reads,"\\|")[[x]][1]))
    df3$i_tumor_ref_reads<-as.numeric(a)
  }
  if (length(which(B!=-1))==0) {
    df3i_tumor_var_reads<-df2$i_tumor_var_reads
  } else {
    b<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_tumor_var_reads,"\\|")[[x]][1]))
    df3$i_tumor_var_reads<-as.numeric(b)
  }
  if (length(which(C!=-1))==0) {
    df3$i_tumor_VAF<-df2$i_tumor_VAF
  } else {
    c<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_tumor_VAF,"\\|")[[x]][1]))
    df3$i_tumor_VAF<-as.numeric(c)
  }
  vaf_dat[[i]]<-data.frame(gene=df$Hugo_Symbol,df3)
  sc[[i]] <- sciClone(vafs=df3, sampleNames=unique(df$Tumor_Sample_Barcode), useSexChrs=FALSE, copyNumberCalls=NULL, copyNumberMargins=0.4, minimumDepth=0, verbose=1, doClusteringAlongMargins=F, plotIntermediateResults = plotIntermediateResults,maximumClusters = MAXCLUSTER,annotation = annot, clusterMethod = 'gaussian.bmm')
  rm(A,B,C)
  rm(df,df1,df2,df3,annot)
  cat(i,"\n")
}

TCGA_patient<-substr(TCGA_barcode,1,12)
names(sc)<-TCGA_patient
names(vaf_dat)<-TCGA_patient
vaf_dat_gene<-lapply(1:length(vaf_dat),function(x) as.character(vaf_dat[[x]]$gene))
all_gene<-unique(unlist(vaf_dat_gene))
pat_vaf<-matrix(0,length(TCGA_patient),length(all_gene))
rownames(pat_vaf)<-TCGA_patient
colnames(pat_vaf)<-all_gene
for (i in 1:length(TCGA_patient)) {pat_vaf[i,vaf_dat_gene[[i]]]<-as.numeric(vaf_dat[[i]][,6])}

m<-NULL;for (i in 1:length(sc)) {m[i]<-length(sc[[i]])}
ig_id<-which(m==0)
sc.1<-sc[-ig_id]
vaf_dat.1<-vaf_dat[-ig_id]

cp_mut_mat<-list()
vaf_dat1_gene<-list()
for (i in 1:length(vaf_dat.1)) {
  A<-vaf_dat.1[[i]]
  p<-which(A$Chromosome=="X" | A$Chromosome=="Y")
  if (length(p)==0) {B<-A} else {B<-A[-p,]}
  vaf_dat1_gene[[i]]<-B$gene
  C<-sc.1[[i]]@clust$cluster.assignments
  M<-matrix(0,nrow(B),length(levels(as.factor(C))))
  rownames(M)<-B$gene
  for (k in 1:nrow(B)) {M[k,C[k]]<-1}
  cp_mut_mat[[i]]<-M
  rm(A,B,C,M)
  cat(i,"\n")
}
names(cp_mut_mat)<-names(sc.1)

path<-"~/Documents/DeepDRpaper-version2/versionHongKong2018/sciClone/SARC.Mutation_Packager_Oncotated_Calls.Level_3"
files <- dir(path, recursive=TRUE, full.names=TRUE, pattern="\\.txt$")

vaf_column<-c("Chromosome","Start_position","t_ref_count","t_alt_count","i_tumor_vaf")
TCGA_barcode<-NULL
library(sciClone)
MAXCLUSTER=10
plotIntermediateResults = 0

sc<-list()
vaf_dat<-list()
for (i in 1:length(files)) {
  df<-read_delim(files[i],  "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)
  TCGA_barcode[i]<-unique(df$Tumor_Sample_Barcode)
  annot = df[,c(1,5,6)]
  df1<-df[,vaf_column]
  df2<-data.frame(df1)
  df3<-df2
  A<-unlist(gregexpr("\\|",df2$t_ref_count))
  B<-unlist(gregexpr("\\|",df2$t_alt_count))
  C<-unlist(gregexpr("\\|",df2$i_tumor_vaf))
  if (length(which(A!=-1))==0) {
    df3$t_ref_count<-df2$t_ref_count
  } else {
    a<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$t_ref_count,"\\|")[[x]][1]))
    df3$t_ref_count<-as.numeric(a)
  }
  if (length(which(B!=-1))==0) {
    df3$t_alt_count<-df2$t_alt_count
  } else {
    b<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$t_alt_count,"\\|")[[x]][1]))
    df3$t_alt_count<-as.numeric(b)
  }
  if (length(which(C!=-1))==0) {
    df3$i_tumor_vaf<-df2$i_tumor_vaf
  } else {
    c<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_tumor_vaf,"\\|")[[x]][1]))
    df3$i_tumor_vaf<-as.numeric(c)
  }
  vaf_dat[[i]]<-data.frame(gene=df$Hugo_Symbol,df3)
  sc[[i]] <- sciClone(vafs=df3, sampleNames=unique(df$Tumor_Sample_Barcode), useSexChrs=FALSE, copyNumberCalls=NULL, copyNumberMargins=0.4, minimumDepth=0, verbose=1, doClusteringAlongMargins=F, plotIntermediateResults = plotIntermediateResults,maximumClusters = MAXCLUSTER,annotation = annot, clusterMethod = 'gaussian.bmm')
  rm(A,B,C)
  rm(df,df1,df2,df3,annot)
  cat(i,"\n")
}

TCGA_patient<-substr(TCGA_barcode,1,12)
names(sc)<-TCGA_patient
names(vaf_dat)<-TCGA_patient
vaf_dat_gene<-lapply(1:length(vaf_dat),function(x) as.character(vaf_dat[[x]]$gene))
all_gene<-unique(unlist(vaf_dat_gene))
pat_vaf<-matrix(0,length(TCGA_patient),length(all_gene))
rownames(pat_vaf)<-TCGA_patient
colnames(pat_vaf)<-all_gene
for (i in 1:length(TCGA_patient)) {pat_vaf[i,vaf_dat_gene[[i]]]<-as.numeric(vaf_dat[[i]][,6])}

m<-NULL;for (i in 1:length(sc)) {m[i]<-length(sc[[i]])}
ig_id<-which(m==0)
sc.1<-sc[-ig_id]
vaf_dat.1<-vaf_dat[-ig_id]

cp_mut_mat<-list()
vaf_dat1_gene<-list()
for (i in 1:length(vaf_dat.1)) {
  A<-vaf_dat.1[[i]]
  p<-which(A$Chromosome=="X" | A$Chromosome=="Y")
  if (length(p)==0) {B<-A} else {B<-A[-p,]}
  vaf_dat1_gene[[i]]<-B$gene
  C<-sc.1[[i]]@clust$cluster.assignments
  M<-matrix(0,nrow(B),length(levels(as.factor(C))))
  rownames(M)<-B$gene
  for (k in 1:nrow(B)) {M[k,C[k]]<-1}
  cp_mut_mat[[i]]<-M
  rm(A,B,C,M)
  cat(i,"\n")
}
names(cp_mut_mat)<-names(sc.1)

path<-"~/Documents/DeepDRpaper-version2/versionHongKong2018/sciClone/TGCT.Mutation_Packager_Oncotated_Calls.Level_3"
files <- dir(path, recursive=TRUE, full.names=TRUE, pattern="\\.txt$")

vaf_column<-c("Chromosome","Start_position","i_NTotCov_SOL","i_TTotCov_SOL","i_TVaf_SOL")
library(sciClone)
MAXCLUSTER=10
plotIntermediateResults = 0

TCGA_barcode<-NULL
sc<-list()
vaf_dat<-list()
for (i in 1:length(files)) {
  df<-read_delim(files[i],  "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)
  TCGA_barcode[i]<-unique(df$Tumor_Sample_Barcode)
  annot = df[,c(1,5,6)]
  df1<-df[,vaf_column]
  df2<-data.frame(df1)
  df3<-df2
  A<-unlist(gregexpr("\\|",df2$i_NTotCov_SOL))
  B<-unlist(gregexpr("\\|",df2$i_TTotCov_SOL))
  C<-unlist(gregexpr("\\|",df2$i_TVaf_SOL))
  if (length(which(A!=-1))==0) {
    df3$i_NTotCov_SOL<-df2$i_NTotCov_SOL
  } else {
    a<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_NTotCov_SOL,"\\|")[[x]][1]))
    df3$i_NTotCov_SOL<-as.numeric(a)
  }
  if (length(which(B!=-1))==0) {
    df3$i_TTotCov_SOL<-df2$i_TTotCov_SOL
  } else {
    b<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_TTotCov_SOL,"\\|")[[x]][1]))
    df3$i_TTotCov_SOL<-as.numeric(b)
  }
  if (length(which(C!=-1))==0) {
    df3$i_TVaf_SOL<-100*df2$i_TVaf_SOL
  } else {
    c<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_TVaf_SOL,"\\|")[[x]][1]))
    df3$i_TVaf_SOL<-100*as.numeric(c)
  }
  vaf_dat[[i]]<-data.frame(gene=df$Hugo_Symbol,df3)
  sc[[i]] <- sciClone(vafs=df3, sampleNames=unique(df$Tumor_Sample_Barcode), useSexChrs=FALSE, copyNumberCalls=NULL, copyNumberMargins=0.4, minimumDepth=0, verbose=1, doClusteringAlongMargins=F, plotIntermediateResults = plotIntermediateResults,maximumClusters = MAXCLUSTER,annotation = annot, clusterMethod = 'gaussian.bmm')
  rm(A,B,C)
  rm(df,df1,df2,df3,annot)
  cat(i,"\n")
}

TCGA_patient<-substr(TCGA_barcode,1,12)
u_TCGApatient<-unique(TCGA_patient)
sc_p<-lapply(1:length(u_TCGApatient),function(x) sc[which(TCGA_patient %in% u_TCGApatient[x])])
names(sc_p)<-u_TCGApatient
vaf_dat_p<-lapply(1:length(u_TCGApatient),function(x) vaf_dat[which(TCGA_patient %in% u_TCGApatient[x])])
names(vaf_dat_p)<-u_TCGApatient
vaf_datp_gene<-lapply(1:length(vaf_dat_p),function(x) as.character(vaf_dat_p[[x]][[1]][,1]))
all_gene<-unique(unlist(vaf_datp_gene))
pat_vaf<-matrix(0,length(u_TCGApatient),length(all_gene))
rownames(pat_vaf)<-u_TCGApatient
colnames(pat_vaf)<-all_gene
for (i in 1:length(u_TCGApatient)) {pat_vaf[i,vaf_datp_gene[[i]]]<-as.numeric(vaf_dat_p[[i]][[1]][,6])}

m<-NULL;for (i in 1:length(sc_p)) {m[i]<-length(sc_p[[i]][[1]])}
ig_id<-which(m==0) # can not get clone combinations because leass of ten samples in vaf dara

cp_mut_mat<-list()
vaf_dat1_gene<-list()
for (i in 1:length(vaf_dat_p)) {
  A<-vaf_dat_p[[i]][[1]]
  p<-which(A$Chromosome=="X" | A$Chromosome=="Y")
  if (length(p)==0) {B<-A} else {B<-A[-p,]}
  vaf_dat1_gene[[i]]<-B$gene
  C<-sc_p[[i]][[1]]@clust$cluster.assignments
  M<-matrix(0,nrow(B),length(levels(as.factor(C))))
  rownames(M)<-B$gene
  for (k in 1:nrow(B)) {M[k,C[k]]<-1}
  cp_mut_mat[[i]]<-M
  rm(A,B,C,M)
  cat(i,"\n")
}
names(cp_mut_mat)<-names(sc_p)

path<-"~/Documents/DeepDRpaper-version2/versionHongKong2018/sciClone/UCEC.Mutation_Packager_Oncotated_Calls.Level_3"
files <- dir(path, recursive=TRUE, full.names=TRUE, pattern="\\.txt$")
vaf_column<-c("Chromosome","Start_position","i_normal_depth","i_tumor_depth","i_tumor_vaf")

library(sciClone)
MAXCLUSTER=10
plotIntermediateResults = 0

TCGA_barcode<-NULL
sc<-list()
vaf_dat<-list()
for (i in 1:length(files)) {
  df<-read_delim(files[i],  "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)
  TCGA_barcode[i]<-unique(df$Tumor_Sample_Barcode)
  annot = df[,c(1,5,6)]
  df1<-df[,vaf_column]
  df2<-data.frame(df1)
  df3<-df2
  A<-unlist(gregexpr("\\|",df2$i_normal_depth))
  B<-unlist(gregexpr("\\|",df2$i_tumor_depth))
  C<-unlist(gregexpr("\\|",df2$i_tumor_vaf))
  if (length(which(A!=-1))==0) {
    df3$i_normal_depth<-as.numeric(df2$i_normal_depth)
  } else {
    a<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_normal_depth,"\\|")[[x]][1]))
    df3$i_normal_depth<-as.numeric(a)
  }
  if (length(which(B!=-1))==0) {
    df3$i_tumor_depth<-as.numeric(df2$i_tumor_depth)
  } else {
    b<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_tumor_depth,"\\|")[[x]][1]))
    df3$i_tumor_depth<-as.numeric(b)
  }
  if (length(which(C!=-1))==0) {
    df3$i_tumor_vaf<-as.numeric(substr(df2$i_tumor_vaf,1,nchar(df2$i_tumor_vaf)-1))
  } else {
    c<-unlist(lapply(1:nrow(df2),function(x) strsplit(df2$i_tumor_vaf,"\\|")[[x]][1]))
    cc<-as.numeric(substr(c,1,nchar(c)-1))
    f<-100*cc[which(cc<1)]
    cc[which(cc<1)]<-f
    df3$i_tumor_vaf<-cc
  }
  vaf_dat[[i]]<-data.frame(gene=df$Hugo_Symbol,df3)
  df4 <- na.omit(df3)
  sc[[i]] <- sciClone(vafs=df4, sampleNames=unique(df$Tumor_Sample_Barcode), useSexChrs=FALSE, copyNumberCalls=NULL, copyNumberMargins=0.4, minimumDepth=0, verbose=1, doClusteringAlongMargins=F, plotIntermediateResults = plotIntermediateResults,maximumClusters = MAXCLUSTER,annotation = annot, clusterMethod = 'gaussian.bmm')
  rm(A,B,C)
  rm(df,df1,df2,df3,annot)
  cat(i,"\n")
}

TCGA_patient<-substr(TCGA_barcode,1,12)
names(sc)<-TCGA_patient
names(vaf_dat)<-TCGA_patient
vaf_dat_gene<-lapply(1:length(vaf_dat),function(x) as.character(vaf_dat[[x]]$gene))
all_gene<-unique(unlist(vaf_dat_gene))
pat_vaf<-matrix(0,length(TCGA_patient),length(all_gene))
rownames(pat_vaf)<-TCGA_patient
colnames(pat_vaf)<-all_gene
for (i in 1:length(TCGA_patient)) {pat_vaf[i,vaf_dat_gene[[i]]]<-as.numeric(vaf_dat[[i]][,6])}

m<-NULL;for (i in 1:length(sc)) {m[i]<-length(sc[[i]])}
ig_id<-which(m==0)

cp_mut_mat<-list()
vaf_dat_gene<-list()
for (i in 1:length(vaf_dat)) {
  A<-vaf_dat[[i]]
  p<-which(A$Chromosome=="X" | A$Chromosome=="Y")
  if (length(p)==0) {B<-A} else {B<-A[-p,]}
  vaf_dat_gene[[i]]<-B$gene
  C<-sc[[i]]@clust$cluster.assignments
  M<-matrix(0,nrow(B),length(levels(as.factor(C))))
  rownames(M)<-B$gene
  for (k in 1:nrow(B)) {M[k,C[k]]<-1}
  cp_mut_mat[[i]]<-M
  rm(A,B,C,M)
  cat(i,"\n")
}
names(cp_mut_mat)<-names(sc)

#generate patient clusters for CESC patients
load("~/Documents/DeepDRpaper-version2/versionHongKong2018/DeepDRc/TCGA-CESC-response-vaf-data.RData")
rm(list = setdiff(ls(),c("genome_wustl_edu_CESC_IlluminaGA_DNASeq_curated_Level_2_1_0_0_somatic")))
sample_ls<-levels(as.factor(genome_wustl_edu_CESC_IlluminaGA_DNASeq_curated_Level_2_1_0_0_somatic$Tumor_Sample_Barcode))
vaf_column<-c("Chromosome","Start_Position","normal_ref_reads","tumor_ref_reads","tumor_vaf" )
MAXCLUSTER=10
plotIntermediateResults = 0
sc<-list()
vaf_dat<-list()
for (i in 1:length(sample_ls)) {
  df1<-genome_wustl_edu_CESC_IlluminaGA_DNASeq_curated_Level_2_1_0_0_somatic[which(genome_wustl_edu_CESC_IlluminaGA_DNASeq_curated_Level_2_1_0_0_somatic$Tumor_Sample_Barcode %in% sample_ls[i]),vaf_column]
  df2<-data.frame(df1)
  annot = genome_wustl_edu_CESC_IlluminaGA_DNASeq_curated_Level_2_1_0_0_somatic[which(genome_wustl_edu_CESC_IlluminaGA_DNASeq_curated_Level_2_1_0_0_somatic$Tumor_Sample_Barcode %in% sample_ls[i]),c(1,5,6)]
  sc[[i]] <- sciClone(vafs=df2, sampleNames=sample_ls[i], useSexChrs=FALSE, copyNumberCalls=NULL, copyNumberMargins=0.4, minimumDepth=0, verbose=1, doClusteringAlongMargins=F, plotIntermediateResults = plotIntermediateResults,maximumClusters = MAXCLUSTER,annotation = annot, clusterMethod = 'gaussian.bmm')
  vaf_dat[[i]]<-data.frame(gene=annot$Hugo_Symbol,df1)
  cat(i,"\n")
}

TCGA_patient<-substr(sample_ls,1,12)
names(sc)<-TCGA_patient
names(vaf_dat)<-TCGA_patient

m<-NULL;for (i in 1:length(sc)) {m[i]<-length(sc[[i]])}
ig_id<-which(m==0)
sc.1<-sc[-ig_id]
vaf.dat1<-vaf_dat[-ig_id]
cp_mut_mat<-list()
vaf_dat1_gene<-list()
for (i in 1:length(vaf.dat1)) {
  A<-vaf.dat1[[i]]
  p<-which(A$Chromosome=="X" | A$Chromosome=="Y")
  if (length(p)==0) {B<-A} else {B<-A[-p,]}
  vaf_dat1_gene[[i]]<-B$gene
  C<-sc.1[[i]]@clust$cluster.assignments
  M<-matrix(0,nrow(B),length(levels(as.factor(C))))
  rownames(M)<-B$gene
  for (k in 1:nrow(B)) {M[k,C[k]]<-1}
  cp_mut_mat[[i]]<-M
  rm(A,B,C,M)
  cat(i,"\n")
}
names(cp_mut_mat)<-names(sc.1)
