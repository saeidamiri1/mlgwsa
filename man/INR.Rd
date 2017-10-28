\name{INR}
\alias{INR}
\title{INfluential Range}
\description{
Finding Influential Ranges that might be the parts of genes or chromosmes.
}

\references{
Amiri, S., Ivo, D. (2017). Information Theoretic Approach for Genome-Wide Association Study of Parkinson's Disease.
}

\seealso{
See Also as \code{\link{HFB}}.
}

\examples{
### Suppose hfbs is obtained via HFB, the plot of
### the difference between the sequential points can
### generated using the following codes
#ts.plot(diff(hfbs),xlab="diff", ylab="value")

#hfbs3<-which(diff(hfbs)<3)+1
#ts.plot(diff(hfbs[hfbs3]),xlab="diff", ylab="value")
}

\keyword{Function}
