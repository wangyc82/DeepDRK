#this procedure is for polycolonal model
#suppose you have got the clones for each patient aready
# the colones for each patient is either from like sciClone software or by your own knowledge
# the mono_score is the monoclonal model prediction score with row as the drug
# drug_tox is the toxxicity for test drug
# cutoff is for distinguish sensitives from resistants
DeepDRc<-function(mono_score,cutoff,drug_tox)
{
  prediction_score<-mono_score*exp(-0.01*drug_tox)
  
  Slab<-which(prediction_score>cutoff,arr.ind = T)
  Rlab<-which(prediction_score<(1-cutoff),arr.ind = T)
  
  F<-matrix(1e+5,nrow(prediction_score),ncol(prediction_score))
  rownames(F)<-rownames(prediction_score)
  colnames(F)<-colnames(prediction_score)
  
  for (i in 1:nrow(Slab)) {F[Slab[i,1],Slab[i,2]]<-1}
  for (i in 1:nrow(Rlab)) {F[Rlab[i,1],Rlab[i,2]]<-0}
  
  Scount<-apply(F,1,function(x) length(which(x==1)))
  ll<-which(Scount==max(Scount))
  
  if (length(ll)>1) {seed<-F[ll[1],];sl<-ll[1]} else {seed<-F[ll,];sl<-ll}
  
  sSl<-which(seed==1)
  sRl<-which(seed==0)
  rl<-seq(1,nrow(F),1)[-sl]
  
  Fs1<-F[rl,sSl]
  Fs1_count<-apply(Fs1,1,function(x) length(which(x==1)))
  Fs1_sl<-which(Fs1_count==max(Fs1_count))
  Fr1<-F[rl,sRl]
  if (length(sRl)==1) {Fr1_rl<-which(Fr1[names(Fs1_sl)]==1)} else {
    Fr1_count<-apply(Fr1,1,function(x) length(which(x==1)))
    Fr1_rl<-which(Fr1_count[names(Fs1_sl)]==max(Fr1_count[names(Fs1_sl)]))
  }
  
  final_DC<-c(rownames(F)[sl],names(Fr1_rl))
  if (length(final_DC)==1) {final_Res<-F[final_DC,]} else {
    final_Res<-colSums(F[final_DC,])
    final_Res[which(final_Res>=1 & final_Res<1e+5)]<-1
  }
  DeepDRc<-final_DC
  #the final response status for patient depends on the cutoffs here
  #if let cutoff=0.7, means if let 70% clones response to that drug combination, then the final decision is response
  #res_status<-ifelse(length(which(final_Res==1))>=0.7*length(final_Res),"responder","non-responder")
}






