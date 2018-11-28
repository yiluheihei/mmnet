\name{updateKEGGPathway}
\alias{updateKEGGPathway}
\title{
Updating and saving the reference (KEGG metabolic pathway) data 
}
\description{
Pasing the  KEGG Markup Language (KGML) file of KEGG metabolic pathway,  obtaining the KO and corresponding metablites and constructing the reference network.
}
\usage{
updateKEGGPathway(path = Sys.getenv("HOME"))
}
\arguments{
  \item{path}{
A character, the path that reference data saved, default is the user's HOME directory.
}
}
\details{
The KEGG metabolic pathway is the best organized part of the KEGG/PATHWAY database, each metabolic pathway can be viewed as a network of KOs. The KEGG Markup Language (KGML) is an exchange format of the KEGG pathway maps, which is converted from internally used KGML+ (KGML+SVG) format. It will retrieves all KOs, products and substrates from the KGML files of KEGG metabolic pathway. There is an initial reference data saved in the subdirectory data of this package. However, KEGG metabolic pathway updated frequently, this function will update the data which contains KO, substrates, products and refence metabolic network. . For more details, see \code{\link{constructMetabolicNetwork}}.

This function updates the KEGG reference data, and saving in the dir\emph{.mmnet} under users Specified, default is the user's home directory. Date, user name and the R VERSION also saved in a log file. 
}

\references{
\url{http://www.kegg.jp/kegg/xml/docs/}
}
\note{
For this function, download KGML file and reference mebabolic constructing may take a long time 
}

\seealso{
\code{\link{constructMetabolicNetwork}},\code{\link{RefDbcache}}
}
\examples{
#updateMetabolicNetwork()
}

