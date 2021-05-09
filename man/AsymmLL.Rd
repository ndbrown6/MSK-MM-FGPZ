\name{AsymmLL}
\alias{AsymmLL}
\title{
Asymmetric Log-Likelihood of mosaic mutations
}
\description{
Asymmetric Log-Likelihood of mosaic mutations
}
\usage{
AsymmLL(m, c, a, nb)
}
\arguments{
\item{m}{Integer vector of number of alternate allele counts}
\item{c}{Integer vector of total number of reference and alternate counts}
\item{a}{Numeric vector of expected VAF multipliers. The length of
	 \eqn{a} should be equal to \eqn{nb}}
\item{nb}{Integer indicating how many branches the model includes.
	  \eqn{nb = 4} or \eqn{nb = 15}}
}
\value{
Invisibly the value of the likelihood
}
\author{
\cr
David Brown\cr
\email{brownd7@mskcc.org}\cr
}
\keyword{}
