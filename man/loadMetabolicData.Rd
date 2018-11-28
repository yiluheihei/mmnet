\name{loadMetabolicData}
\alias{loadMetabolicData}
\title{
Metabolic data loading
}
\description{
This function load the updated reference metabolic data into workspace.
}
\usage{
loadMetabolicData(path = Sys.getenv("HOME"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
 character, the path that reference data saved, default is the user's HOME directory. If does not exits, prestored reference metabolic data (RefDbCache) in \emph{mmnet} will be loaded.
}
}
\seealso{\code{link{updateKEGGPathway}}, \code{\link{constructMetabolicNetwork}}
}
\examples{
##
# loadMetabolicData(path = Sys.getenv("HOME"))
}