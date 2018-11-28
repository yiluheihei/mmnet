\name{constructMetabolicNetwork}
\alias{constructMetabolicNetwork}
\title{
Metabolic network construction of microbiome
}
\description{
Metabolic network construction according to the metabolites among enzymes
}
\usage{
constructMetabolicNetwork(path = Sys.getenv("HOME"))
}
\arguments{
 \item{path}{
 character, the path that reference data saved, default is the user's HOME directory. If does not exits, prestored reference metabolic data (RefDbCache) in \emph{mmnet} will be loaded.
}
}

\value{
referce metabolic network of igraph
}

\seealso{\code{\link{loadMetabolicData}}}

\examples{
## may take a long time
# data(RefDbcache)
# constructMetabolicNetwork(RefDbcache)
}
