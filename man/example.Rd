\name{example}
\alias{example}
\docType{data}
\title{This data set includes the calculation run on the Parkinson data.
}
\description{
}
\usage{data("example")}
\format{
  A dataset includes some calculation done on Parkinson data.
  \describe{
    \item{\code{PDfeexample}}{result of running IBC on SNPdatafemPD}
    \item{\code{COfeexample}}{result of running IBC on SNPdatafemCO}
    \item{\code{dt1example}}{result of running Distwholeon informative variables}
    \item{\code{vsfeentRexample}}{result of HBF}
    \item{\code{dt2example}}{result of running Distwhole on HBF}
  }
}
\details{
}
\source{
}
\references{
}
\examples{
data("example")
### Assume SNPdatafemCO and SNPdatafemPD are the SNPs of Parkinson- free and
### Parkinson disease females, respectively.
# data(SNPdatafemPD)
# data(SNPdatafemCO)
### We are not allowed to publish the data publically, so we created
### synthetic data SNPgroupC and SNPgroupD.

### Calculate the uncertainty of variables
#PDfeexample<-IBC(SNPdatafemPD)
#COfeexample<-IBC(SNPdatafemCO)

### Find the informative variables
vsfeent1example<-which(COfeexample>PDfeexample)
length(vsfeent1example)

### Calculate the distance between the data on the informative variables
#dt1example<-Distwhole(SNPdatafemPD,SNPdatafemCO,vsfeent1example)
### Calculate hierarchical dendrogram
hcdt1example<-hclust(as.dist(dt1example), method = "average", members = NULL)

### Define the label of data
##labfem<-c(rep("P",dim(SNPdatafemPD)[1]),rep("C",dim(SNPdatafemCO)[1]))
labfem<-c(rep("P",40),rep("C",101))

### plot the dendrogram
plot(hcdt1example,label=labfem)

Alpha<-0.80
hfbexample<-HFB(vsfeentRexample,Alpha)
length(hfbexample)

#### the folowing codes gives the dengrogram on the informative variables
#### obtained via HFB
#dt2example<-Distwhole(SNPdatafemPD,SNPdatafemCO,hfbexample)

### Calculate hierarchial dendrogram
hcdt2<-hclust(as.dist(dt2example), method = "average", members = NULL)
### plot the dendrogram
plot(hcdt2,label=labfem)

ts.plot(diff(hfbexample),xlab="diff", ylab="value")
hfbs3<-which(diff(hfbexample)<3)+1
ts.plot(diff(hfbexample[hfbs3]),xlab="diff", ylab="value")
}
\keyword{datasets}
