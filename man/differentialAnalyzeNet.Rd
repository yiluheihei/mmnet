\name{differentialAnalyzeNet }
\alias{differentialAnalyzeNet }

\title{
the KO difference abundance of specific state
}
\description{
The function computes the difference abundance with three methods and show the differential abundance comparison metabolic networks.
}
\usage{
differentialAnalyzeNet (ssns, sample.state, method = c("OR", "rank", "JSD"), cutoff,  p.value = TRUE, Visualization = TRUE)
}
\arguments{
  \item{ssns}{
the state specific network for differential comparison 
}
  \item{sample.state}{
a character vector, represents the state of sample. Apparently, its length is equal to  number of samples 
  }
  \item{method}{
method used to mesure the difference abundance, more description in details
}
\item{cutoff}{
numeric vector with length 2, specifies the difference thereshold
}
  \item{p.value}{
logical, whether p value is caculated, Wilcoxon rank-sum test to limit enzymes that enriched or depleted
}
  \item{Visualization}{
logical, whether plot the differential metabolic network
}
}
\details{
This function only applicable to compare the abundance profile between samples in different states (e.g., lean vs. obese). States is added to the abundance profile 
as the name attribute. Abundance were normalized within each sample to represent the relative abundance of each enzyme (KO) in each sample before comparison, thereby accounting for differences in sampling depth.


The selection of method to measure the enrichment is a non-trival component of comparative metagenomic analysis, since statistical results can be heavily influenced by missing data, sample size, data magnitute, normalization technique, or assumptions of the distribution of values. 

Odds ratio (OR) is a  widly recognized measurement as well as has a number of benefits. The differential abundance score of each enzyme, defined as abs[log2(OR)], provides a measure of the extent to which an enzyme's abundance differs in samples from a given state, relative to samples in another state. 

RANK: first ranking the enzymes within each sample from most abundant to least abundant, then measure the difference between the mean reank of each enzyme between samples in different states.

 JSD is the Jensen-Shannon divergence measure to quantify the similarity bewteen the distributions of enzyme in samples in different state. Note it only suitable for a large enough
 sample size.

cutoff: more description see \code{\link{showMetagenomicNet}}
}
\value{
a igraph with the diference abundance as a node attribute
}