\name{HFB}
\alias{HFB}
\title{High Frequency via Bootstrap
}
\description{
}
\usage{
HFB(x,alpha)
}
\arguments{
  \item{x}{matrix including the bootstrap of influence variables}
  \item{alpha}{minimum relative frquency tof influence variables.}
}
\details{
To have more reliable informative variables, we implement the HFB that run the proposed method on the resample from data.
}

\references{
Amiri, S., Ivo, D. (2017). Information Theoretic Approach for Genome-Wide Association Study of Parkinson's Disease.
}

\seealso{
}
\examples{
### obtaine the size of dataset.
#nD<-dim(SNPgroupD)[1]
#nC<-dim(SNPgroupC)[1]

### number of bootstrap
#B<-500
#vsfeentR<-list()
### the following loop achieves the bootsrap and obtain the informative variables
#for(b in 1:B){

#  sxx1<-sort(unique(sample(nD,nD,replace=TRUE)))
#  sxx2<-sort(unique(sample(nC,nC,replace=TRUE)))

#  SNPgroupDR<-SNPgroupD[ sxx1,]
#  SNPgroupCR<-SNPgroupC[ sxx2,]
#  DfeR<-IBC(SNPgroupDR)
#  CfeR<-IBC(SNPgroupCR)

#  vsfeentR[[b]]<-which(CfeR/2>DfeR)
#}

#Alpha<-0.8
#hfb<-HFB(vsfeentR,Alpha)
#length(hfb)

#### the following codes give the dengrogram on the informative variables
#### obtained via HFB
#dt1<-Distwhole(SNPgroupD,SNPgroupC,vsfeent1)
### Calculate hierarchial dendrogram
#hcdt1<-hclust(as.dist(dt1), method = "average", members = NULL)
### plot the dendrogram
#plot(hcdt1,label=labfem)
}

\keyword{Function}
