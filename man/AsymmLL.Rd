\name{AsymmLL}
\alias{AsymmLL}
\title{
Asymmetric Log-Likelihood of mosaic mutations
}
\description{
Computes the Log-Likelihood of mosaic mutations given asymmetric
doubling at given cell generations
}
\usage{
AsymmLL(m, c, a, nb)
}
\arguments{
\item{m}{Integer vector of alternate allele counts}
\item{c}{Integer vector of total reference and alternate counts}
\item{a}{Numeric vector of expected VAF multipliers}
\item{nb}{Integer indicating how many cell generations the model includes}
}
\value{
Invisibly a list with the value of the Log-Likelihood
}
\author{
\cr
David Brown\cr
\email{brownd7@mskcc.org}\cr
}
\keyword{}
