\name{HFB}
\alias{HFB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ High Frequency via Bootstrap
%%  ~~function to do ... ~~
}
\description{

}
\usage{
HFB(x,alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{matrix including the bootstrap of influence variables}
  \item{alpha}{minimum relative frquency tof influence variables.}
}
\details{
To have more reliable informative variable, we implement the HFB that run the proposed method on the resample from data and take the more often informative variables.
}

\references{
Amiri, S., Ivo, D. (2015). Genome Wide-Association Study with application of Parkinson's Disease.
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
### obtaine the size of dataset.
#nPD<-dim(SNPdatafemPD)[1]
#nCO<-dim(SNPdatafemCO)[1]

### number of bootstrap
#B<-500
#vsfeentR<-list()
### the following loop achieves the bootsrap and obtain the informative variables
#for(b in 1:B){

#  sxx1<-sort(unique(sample(nPD,nPD,replace=TRUE)))
#  sxx2<-sort(unique(sample(nCO,nCO,replace=TRUE)))

#  SNPdatafemPDR<-SNPdatafemPD[ sxx1,]
#  SNPdatafemCOR<-SNPdatafemCO[ sxx2,]
#  PDfeR<-IBC(SNPdatafemPDR)
#  COfeR<-IBC(SNPdatafemCOR)

#  vsfeentR[[b]]<-which(COfeR/2>PDfeR)
#}

#Alpha<-0.8
#hfb<-HFB(vsfeentR,Alpha)
#length(hfb)

#### the following codes give the dengrogram on the informative variables
#### obtained via HFB
#dt1<-Distwhole(SNPdatafemPD,SNPdatafemCO,vsfeent1)
### Calculate hierarchial dendrogram
#hcdt1<-hclust(as.dist(dt1), method = "average", members = NULL)
### plot the dendrogram
#plot(hcdt1,label=labfem)
}

\keyword{Function}
