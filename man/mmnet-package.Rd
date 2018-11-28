\name{mmnet-package}
\alias{mmnet-package}
\title{
A metagenomic pipeline for systems biology
}
\description{
The \pkg{mmnet} add-on is an implementation of mapping metagenomic data to the KEGG global metabolic pathway and constructing a systems-level network. 
It is of great help to analysis the various functional properties, including regulation and metabolic functionality of the metagenome.
} 
\details{
This package contains functions for mapping metagenomic data to the KEGG global metabolic pathway and constructing systems-level network. First,  the information of enzymes, corresponding substrates and products is parsed from the KEGG PATYWAY database with function \code{\link{updateKEGGPathway}}. Second, sequence functional ontology annotation with KEGG is preprocessing on MGRAST, KO abundance can be obtained by parse the MGRAST ontology annotation files with function \code{\link{estimateAbundance}}. Finally, global metabolic pathway and systems-level network are made to analyze the metabolic functionality of metagenome with function \code{\link{constructMetabolicNetwork}}, visualization with function \code{\link{showMetagenomicNet}}
}
