\name{listMgrastProject}
\alias{listMgrastProject}
\title{
List all the private projects you summitted to MGRAST
}
\description{
List all your private projects information on MGRAST, the information included job id, MGRAST project id and  MGRAST id ("Accession numbers") which completed calculation. 
}
\usage{
listMgrastProject(login.info)
}
\arguments{
  \item{login.info}{
a list of login info generated after user login into MGRAST which consists of webkey of the file, the invlid time of your webkey, cookie, websession and curlhandle. It is unique for each users
}
}
\details{
Time MGRAST annotated a metagenome depend on the size of your dataset and and the current server load. In practice the time taken will range between a few hours and a week. 

It will list all the MGRAST metagenome id for each user's project which was completely annotated, it also means the in progress metagenome of a project will not show until annotation complete. Only this metagenome list by this function can user to view the annotation relusts. 
}
\value{
Return all the private projects, a list contains the job name users created while job submission, and MGRAST id ("Accession numbers of metagenome"), MGRAST id will be NULL while sample annotation is not completed, and project id MGRAST assigned corresponding to the job id.
                metagenome.id), "?verbosity=metadata")
}

\examples{
# listMgrastProject(login.info)
}
