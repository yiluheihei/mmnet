\name{preprocessKOMetabolites}
\alias{preprocessKOMetabolites}
\title{
Preprocessing of KO's (enzyme's) metabolites
}
\description{
This function preprocess the KO information to delete the redundant metabolites 
}
\usage{
preprocessKOMetabolites(ko.info)
}
\arguments{
  \item{ko.info}{
KO information, can be obtained by \code{\link{getKOPathwayInfo}}
}
}
\details{
When we get the KO information first, there are reversible reaction, duplicate KOs and metabolites. This function merge all the KO information together and delete the redundant information.
}
\value{
list, is composed of KO, substrates and products
}

\seealso{\code{\link{getKOPathwayInfo}}, \code{\link{updateKEGGPathway}}
}
