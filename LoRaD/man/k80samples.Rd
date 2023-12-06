\name{k80samples}
\docType{data}
\alias{k80samples}
\title{Posterior Samples From k80 Substitution Model}
\description{
	This data set contains 10000 samples of the edge length and transition-transversion rate ratio (kappa) along with the unnormalized log posterior (log.kernel) for data simulated with true edge length 0.2 and true kappa 5.0 and sequence length 200 base pairs.
}
\usage{k80samples}
\format{A data frame with 4 columns (iter, log.kernel, edgelen, kappa) and 10000 rows excluding the header.}
\source{https://doi.org/10.1093/sysbio/syad007}
\references{
  Wang et al. (2023) \emph{LoRaD: marginal likelihood estimation with haste (but no waste)}.
  Systematic Biology 72(3):639-648.
}
\keyword{datasets}