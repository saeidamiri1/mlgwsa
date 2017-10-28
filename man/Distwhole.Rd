\name{Distwhole}
\alias{Distwhole}
\title{Distwhole}
\usage{
Distwhole(x,y,inf)
}
\arguments{
  \item{x}{matrix includes the first dataset.}
  \item{y}{matrix includes the second dataset.}
  \item{inf}{vector of the informative SNP.}
}
\details{
This function calculates the distance between the data sets using the informative variables. First calculate the entropy of variable, the uncertainty of variables, using "IBC", then compute the distance on the the informative variables.
}
\examples{
### We use two synthetic data to explain the codes
#data(SNPgroupC)
#data(SNPgroupD)
### Define the label of data
#labfem<-c(rep("C",dim(SNPgroupC)[1]),rep("D",dim(SNPgroupD)[1]))

### Calculate the uncertainty of variables
#Dfe<-IBC(SNPgroupD)
#Cfe<-IBC(SNPgroupC)

### Find the informaive variables
#vsfeent1<-which(Cfe>Dfe)
#length(vsfeent1)

### Calculate the ditance between the data on the informative variables
#dt1<-Distwhole(SNPgroupD,SNPgroupC,vsfeent1)
### Calculate the hierarchial dendrogram
#hcdt1<-hclust(as.dist(dt1), method = "average", members = NULL)
### plot the dendrogram
#plot(hcdt1,label=labfem)

### The proposed methods can be done with more constraint, c=2
#vsfeent2<-which(Cfe/2>Dfe)
#length(vsfeent2)
#dt2<-Distwhole(SNPgroupD,SNPgroupC,vsfeent2)
#hcdt2<-hclust(as.dist(dt2), method = "average", members = NULL)
#plot(hcdt2,label=labfem)
}
