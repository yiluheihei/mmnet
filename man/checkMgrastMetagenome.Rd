\name{checkMgrastMetagenome}
\alias{checkMgrastMetagenome}
\title{
Check whether the MGRAST annotation is completed
}
\description{
Check the whether MGRAST annotation is comleted of the metagenomic sequence. 
}
\usage{
checkMgrastMetagenome(login.info, metagenome.id, public = TRUE)
}
\arguments{
  \item{login.info}{
login.info a list of login info generated after user login into MGRAST
which consists of webkey of the file, the invlid time of your webkey, cookie,
websession and curlhandle. It is unique for each users. Only needed while checking
the private project.
}
\item{metagenome.id}{
a chracter, metagenome.id the metagenome id MGRAST assigned while user submit 
metagenome sequence to MGRAST.
}
\item{public}{
logical, whether the metagenome of your MGRAST project is public
}
}
\value{
logical, whether the metagenome annotation is completed, it will be NULL while metagenome annotation is not completed.
}
\examples{
checkMgrastMetagenome(metagenome.id = "4440616.3")
}
