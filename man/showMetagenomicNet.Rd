\name{showMetagenomicNet}
\alias{showMetagenomicNet}
\title{
Visualiation of metabolic network 
}
\description{
This function able to plot metabolic network to any R device. 
}
\usage{
showMetagenomicNet(net, mode = c("ref", "ssn", "compared"), method = c("OR","rank", "JSD"), cutoff, ...)
}
\arguments{
  \item{net}{
the graph to visualization, reference metabolic network or sub-network of different samples
}
  \item{mode}{
a string, represents the property of your network, for more description in detalils
}
\item{method}{
method used to mesure the difference abundance, more description in \code{\link{differentialAnalyzeNet}}
}
 \item{cutoff}{
numeric vector with length 2, specifies the difference thereshold, more description in details
}
  \item{\dots}{   
additional plotting parameters. See \code{\link{igraph.plotting}} for the complete list.
}
}
\details{
There are three different metabolic networks: ref - reference global network, ssn - state specific, compared - compared network between different state network. 

More details on  method to measure the difference abundance see \code{\link{differentialAnalyzeNet}}

cutoff: two-fold threshold commonly used in the odds-ratio test, cutoff = c(0.5, 2); enzymes with the 10\% highest calculated differences and the 10\% lowest differences, cutoff = c(0.1, 0.9)
}
\value{
Returns NULL, invisibly.
}

\examples{
# reference network
data(RefDbcache)
showMetagenomicNet(RefDbcache$network,mode="ref")
}

