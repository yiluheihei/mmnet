\name{constructSSN}
\alias{constructSSN}
\title{
construct the state specific network
}
\description{
Given a set of KOs of a sample with a specific state and the reference metabolic network, this function creates and returns sub-network with only the supplied nodes (KOs) and any edges between them.
}
\usage{
constructSSN(abundance)
}
\arguments{
\item{abundance}{ 
abundance of KOs of a given specific state, may a numeric vector which can be obtained with \code{\link{estimateAbundance}} 

a given Biological Observation Matrix (BIOM) format, the ncol of BIOM file represents number of samples. The BIOM file format (canonically pronounced biome) is designed to be a general-use format for representing biological sample by observation contingency tables. More details on BIOM file can be found in \url{http://biom-format.org/}  R package \emph{biom}. There is also a example BIOM file in directory "inst/extdata" of \emph{mmnet}. 

It could be obtained with \code{\link{estimateAbundance}}. Moreover, while samples has been annotated with other tools (e.g. blast to KEGG locally), users could create their own BIOM file represents enzymatic genes abundance.    
}
}
\details{
Reference metabolic network has been saved in the RefDbcache. There are varies sets of KOs and its abundance between metagenomic samples. With this function, users can create a specified sub-network of the reference network for each sample.
}
\value{
new graph object will be return, and the abundance is added as a vertex meta-attributes. The number of graph is equal to the samples. 
}

\seealso{
\code{\link{RefDbcache}}, \code{\link{estimateAbundance}}
}
\examples{
# load reference metabolic network 
data(RefDbcache)

data(anno)
abund <- estimateAbundance(anno[[1]])
subnet  <- constructSSN(abund)

}

