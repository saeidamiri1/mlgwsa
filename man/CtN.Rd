\name{CtN}
\alias{CtN}
\title{Characteristic to Numerical
}
\description{This function converts genetic data (nucleotides) to numeric data
}
\usage{CtN(sdata)}
\arguments{
  \item{sdata}{A SNP data with \{A,T,C,G\}.}
}

\details{
R is faster on the numerical data and storage of data via numerical values takes less capacity. Genomics data includes nucleotides \{A,T,C,G\} and are often saved in Fasta format. This function converts genomics data (nucleotides) to numeric data; converts \{A,T,C,G\} to \{1, 2, 3, 4\} and combine them together.
}

\examples{
# SNPexample
# CtN(SNPexample)
}
\keyword{datasets}
