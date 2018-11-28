\name{saveMetabolicData}
\alias{saveMetabolicData}
\title{
Save a new version of the reference metabolic data
}
\description{
This function saves the updated reference KEGG metabolic data and global network in your computer.
}
\usage{
saveMetabolicData(RefDbcache, path = Sys.getenv("HOME"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RefDbcache}{
a list, the updated KEGG Pathway information, more details in \code{link{preprocessKOMetabolites}} 
}
  \item{path}{
character, the path that reference data saved, default is the user's HOME directory.
}
}
\seealso{
\code{link{preprocessKOMetabolites}} 
}
\examples{
## 
# saveMetabolicData(RefDbcache, path = Sys.getenv("HOME"))
}