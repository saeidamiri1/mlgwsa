\name{INR}
\alias{INR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ INfluential Range
%%  ~~function to do ... ~~
}
\description{
Finding Influential Ranges that migh be parts of genes or chromosmes.
The easiet way is to find the sequential point of indices that is
obtained via HFB.
}

\references{
Amiri, S., Ivo, D. (2015). Genome Wide-Association Study with application of Parkinson's Disease.
}

\seealso{
See Also as \code{\link{HFB}}.
}
\examples{
### Suppose hfbs is obtaine via HFB, the distance
### between the sequential point and plot of them can be obtained via
#ts.plot(diff(hfbs),xlab="diff", ylab="value")

#hfbs3<-which(diff(hfbs)<3)+1
#ts.plot(diff(hfbs[hfbs3]),xlab="diff", ylab="value")
}

\keyword{ Function }
