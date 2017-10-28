CtN<-function(x){
  CN4<-function(XX0){
    XX1<-NULL
    XX1[XX0=='0'| XX0=='I'| XX0=='D' ]<-NA
    XX1[XX0=='A']<-1
    XX1[XX0=='T']<-2
    XX1[XX0=='C']<-3
    XX1[XX0=='G']<-4
    return(XX1)
  }
  x<-t(x)
  data0<-matrix(,nrow=(dim(x)[1]), ncol=(dim(x)[2]))
  for(i in 1:(dim(x)[2])){
    data0[,i]<-CN4(x[,i])
  }
  t(data0)
}
################
################
################
################
################
################
################
################
################
################

ENTO<-function(cc) {
  p0<-table(cc)/sum(table(cc))
  -mean(p0*log(p0))
}


########################
hammingD<-DistM<-function(Xdat){
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


##################
##################




