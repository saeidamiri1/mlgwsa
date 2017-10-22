ENTO<-function(cc) {
  p0<-table(cc)/sum(table(cc))
  -mean(p0*log(p0))
}


########################
DistM<-function(Xdat){
  Len<-dim(Xdat)
  ss<-Len[1]
  dismat<- matrix(1,ncol=ss,ss)
  for(i in 1:ss){
    if (i==ss) break
    for(j in (i:ss)){
      dismat[i,j]<-mean(Xdat[i,]==Xdat[j,],na.rm=T)
    }
  }
  return(1-t(dismat))
}


IBC<-function(data0){

  homogout<-apply(data0,2,ENTO)

return(homogout)

}



Distwhole<-function(SNPfemale1,SNPfemale2,vs1feent){
  SNPfemale<-rbind(SNPfemale1,SNPfemale2)
  outdist<-DistM(SNPfemale[,vs1feent])
  return(outdist)
}



#######################

HFB<-function(x,alpha){
  wsw<-unlist(x)
  B<-length(x)*alpha
  #table(wsw)
  ux <- unique(wsw)
  return(sort(ux[which(tabulate(match(wsw, ux))>B)]))
}
