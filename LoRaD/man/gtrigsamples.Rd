\name{gtrigsamples}
\docType{data}
\alias{gtrigsamples}
\title{Posterior Samples From GTR+I+G Substitution Model}
\description{
	This data set contains 10000 samples of parameters from a fixed tree topology and using the GTR+I+G model: shape parameter of the rate heterogeneity Gamma distribution (alpha), proportion of invariable sites (pinvar), tree length (tree_length), edge length proportions (edge_length_proportions), and nucleotide relative frequencies (pi).
	It also contains the log posterior (Posterior), likelihood (Likelihood), and prior (Prior) as well as functions of parameters, including edge lengths (edgelens) and site rates (site_rates).
	The data are 1296 base pairs long.
	}
\usage{gtrigsamples}
\format{A data frame with 35 columns (Iteration, Posterior, Likelihood, Prior, alpha, edge_length_proportions[1], edge_length_proportions[2], edge_length_proportions[3], edge_length_proportions[4], edge_length_proportions[5], edge_length_proportions[6], edge_length_proportions[7], edgelens[1], edgelens[2], edgelens[3], edgelens[4], edgelens[5], edgelens[6], edgelens[7], er[1], er[2], er[3], er[4], er[5], er[6], pi[1], pi[2], pi[3], pi[4], pinvar, site_rates[1], site_rates[2], site_rates[3], site_rates[4], tree_length) and 10000 rows excluding the header.}
\source{https://doi.org/10.1093/sysbio/syad007}
\references{
  Wang et al. (2023) \emph{LoRaD: marginal likelihood estimation with haste (but no waste)}.
  Systematic Biology 72(3):639-648.
}
\keyword{datasets}
