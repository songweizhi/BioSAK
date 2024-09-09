
spec.gen<-function(comm.tab,niche.width.method="levins",perm.method="quasiswap",n=1000,probs=c(0.025,0.975)){
  require(spaa)
  require(vegan)
  occurrence<-function(x){apply(ceiling(x/max(x)),2,sum)}
  n<-n
  if (niche.width.method=="occurrence") levin.index.real<-occurrence(comm.tab) else levin.index.real<-as.numeric(niche.width(comm.tab,method=niche.width.method))
  names(levin.index.real)<-colnames(comm.tab)
  
  levin.index.simul<-matrix(NA,ncol=dim(comm.tab)[2],nrow=n)
  for (i in 1:n){
    if (niche.width.method=="occurrence") levin.index.simul[i,]<-occurrence(permatswap(comm.tab,perm.method,times=1)$perm[[1]]) else levin.index.simul[i,]<-as.numeric(niche.width(permatswap(comm.tab,perm.method,times=1)$perm,method=niche.width.method))
  }
  colnames(levin.index.simul)<-colnames(comm.tab)
  levin.index.simul<-as.data.frame(levin.index.simul)
  media<-apply(levin.index.simul,2,mean)
  ci<-apply(levin.index.simul,2,quantile,probs=probs)
  resultats<-data.frame(observed=levin.index.real,mean.simulated=media,lowCI=ci[1,],uppCI=ci[2,],sign=NA)
  for (j in 1:dim(resultats)[1]){
    if (resultats$observed[j]>resultats$uppCI[j]) resultats$sign[j]<-"GENERALIST"
    if (resultats$observed[j]<resultats$lowCI[j]) resultats$sign[j]<-"SPECIALIST"
    if (resultats$observed[j]>=resultats$lowCI[j] & resultats$observed[j]<=resultats$uppCI[j]) resultats$sign[j]<-"NON SIGNIFICANT"
  }
  resultats$sign<-as.factor(resultats$sign)
  resultats}

