\name{topologicalAnalyzeNet}
\alias{topologicalAnalyzeNet}
\title{
Analyze the specific state network 
}
\description{
Calculate Topological features of each node and linked with abundance attributes.
}
\usage{
topologicalAnalyzeNet(g, Scatterplot = TRUE, .properties = c("betweennessCentrality", "degree", "clusteringCoefficient", "pageRank"), mode = c("all", "in", "out"), ...)
}
\arguments{
  \item{g}{
a igraph network of a specific state, it can be obtained with function \code{\link{constructSSN}} 
}
  \item{Scatterplot}{
logical, if show the scatterplot between topolagical features and abundance
}
  \item{.properties}{
instead of specifying the topological features in free form via the ... argument, one can specify them as a vector
}
  \item{mode}{
Character string, "out" for out-degree, "in" for in-degree and "all" for the sum of the two.
}
  \item{\dots}{
topological features that user want to calculate.
}
}
\details{
Five most presented vertex topological feartures is provided by this function.

betweenness centrality:  defined by the number of geodesics (shortest paths) going through a vertex , more details in\code{\link{betweenness}} of package igraph

degree: the most basic structural property, represents the number of its adjacent edges.

clustringCofficient: measures the probability that the adjacent vertices of a vertex are connected

pageRank: Calculates the Google PageRank for the vertices in network
}

\value{
a graph that add all the calculated topological features as the meta-attributes
}
\seealso{
\code{\link{constructSSN}}
}
\examples{
##
data(anno)
data(RefDbcache)
KOAbund <- estimateAbundance(anno[[1]])

subnet  <- constructSSN(KOAbund)
topologicalAnalyzeNet(subnet)
}
