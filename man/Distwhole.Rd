\name{Distwhole}
\alias{Distwhole}
\title{Distwhole}
\usage{
Distwhole(x)
}
\arguments{
  \item{x}{matrix includes the uncertainty of data obtained via the proposed method.
%%     ~~Describe \code{x} here~~
}
}
\details{
This function calculates the ditance between the data on the informative variables. First the entropy of variable should be found via the function "IBC", the following example includes complete details.
}
\examples{
### the following library includes the Parkinson's data
#library("ParkinsonData")
#data(SNPdatafemCO)
#data(SNPdatafemPD)
### define the label of data
#labfem<-c(rep("P",dim(SNPdatafemPD)[1]),rep("C",dim(SNPdatafemCO)[1]))

### Calculate the uncertainty of variables
#PDfe<-IBC(SNPdatafemPD)
#COfe<-IBC(SNPdatafemCO)

### Find the informaive variables
#vsfeent1<-which(COfe>PDfe)
#length(vsfeent1)

### Calculate the ditance between the data on the informative variables
#dt1<-Distwhole(SNPdatafemPD,SNPdatafemCO,vsfeent1)
### Calculate hierarchial dendrogram
#hcdt1<-hclust(as.dist(dt1), method = "average", members = NULL)
### plot the dendrogram
#plot(hcdt1,label=labfem)


### The proposed methods can be done with more constraint, c=2
#vsfeent2<-which(COfe/2>PDfe)
#length(vsfeent2)
#dt2<-Distwhole(SNPdatafemPD,SNPdatafemCO,vsfeent2)
#hcdt2<-hclust(as.dist(dt2), method = "average", members = NULL)
#plot(hcdt2,label=labfem)


### the distance might be small, it can be standardized by
#hcdt2b<-hclust(as.dist(dt2/max(dt2)), method = "average", members = NULL)
#plot(hcdt2b,label=labfem)

}
