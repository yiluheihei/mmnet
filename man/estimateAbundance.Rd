\name{estimateAbundance}
\alias{estimateAbundance}
\title{
Getting KO (KEGG ontology) and abundance information
}
\description{
This function parses MGRAST functional annotation files to retrieve the KO of the KEGG, and estimate the abundance of each KO. 
}
\usage{
estimateAbundance(KOAnno)
}
\arguments{
  \item{KOAnno}{
for multiple samples: a list, each element represents MGRAST functional annotation file suffixed with "expand.ontology" of one sample, includes md5 value of the database sequence hit followed by sequence or cluster ID, similarity information, and annotation information. Users can download it in the MGRAST website or with function \code{\link{getMgrastAnnotation}}.

for a single sample: a data frame, represents the MGRAST functional annotation file suffixed with "expand.ontology" of the sample
}
}
\details{
With the annotation file, some values are encoded as index numbers in the annotation files of the argument.
Column fields are as follows:

    1. hit id (md5sum) index number

    2. query id (either fasta ID or cluster ID)

    3. percentage identity

    4. alignment length

    5. e-value

    6. function annotation index number

    7. functional category (ontology) index number

    8. DB source index number

This function is just to get the KEGG ontology from the ontology annotation files. Users can substract reference networks by mapping this annotation to the selected reference database (KEGG PATHWAY).
}
\value{
Return a a given Biological Observation Matrix (BIOM) format, the ncol of the BIOM file is equal to the number of sample, and the colnames of the contingency table represents healthy state of samples. The BIOM file format (canonically pronounced biome) is designed to be a general-use format for representing biological sample by observation contingency tables. More details on BIOM file can be found in \url{http://biom-format.org/}  R package \emph{biom}. There is also a example BIOM file in directory "inst/extdata" of \emph{mmnet}. , the name of the vector is the KO number and elements represents abundance.
}

\seealso{
\code{\link{constructMetabolicNetwork}} 
}

\examples{
#load the test ontology data
data(anno)
KOAbund <- estimateAbundance(anno[[1]])
}

