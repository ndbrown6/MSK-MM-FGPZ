\name{SymmLL}
\alias{SymmLL}
\title{
Symmetric Log-Likelihood of mosaic mutations
}
\description{
Computes the Log-Likelihood of mosaic mutations given symmetric
doubling at all cell generations
}
\usage{
AsymmLL(m, c, nb)
}
\arguments{
\item{m}{Integer vector of alternate allele counts}
\item{c}{Integer vector of total reference and alternate counts}
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