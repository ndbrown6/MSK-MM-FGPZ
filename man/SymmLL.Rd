\name{SymmLL}
\alias{SymmLL}
\title{
Likelihood of vb given symmetric cell divisions
}
\description{
Returns the likelihood of a mosaic mutation occuring on a given model lineage tree vb
}
\usage{
SymmLL(m_j, c_j, b = c(1, 2, 3, 4), method = c("binomial", "betabinomial"), log = c(TRUE, FALSE))
}
\arguments{
\item{m_j}{Integer for the number of alternate allele counts}
\item{c_j}{Interger for the total number of reference and alternate counts}
\item{b}{Integer for the level of the level of the model lineage tree vb}
\item{method}{String to choose between a Binomial and a Beta-Binomial distribution}
\item{log}{Logical to choose between raw or log-transformed probabilities}
}
\value{
Invisibly the value of the likelihood vb
}
\author{
\cr
David Brown\cr
\email{brownd7@mskcc.org}\cr
}
\keyword{}
