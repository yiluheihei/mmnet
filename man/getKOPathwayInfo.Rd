\name{getKOPathwayInfo}
\alias{getKOPathwayInfo}
\title{
Get the KEGG metabolic pathway information
}
\description{
This function phases the KEGG metabolic pathway information, and obtains the KO or enayzme ,reaction, metabolites, for each pathway
}
\usage{
getKOPathwayInfo(koPathway)
}
\arguments{
  \item{koPathway}{
the metabolic KO pathway name, can be obtained with the \emph{KEGGREST} package 
}
}

\value{
list of the KEGG metabolic pathway information
}

\seealso{
\code{\link{updateKEGGPathway}}
}

