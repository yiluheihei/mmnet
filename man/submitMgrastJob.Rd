\name{submitMgrastJob}
\alias{submitMgrastJob}
\title{
Create a job on MGRAST
}
\description{
Submit your file to create a job after data uploaded.   
}
\usage{
submitMgrastJob(login.info,seqfile,new_project)
}
\arguments{
  \item{login.info}{
a list of login info generated after user login into MGRAST which consists of webkey of the file, the invlid time of your webkey, cookie, websession and curlhandle. It is unique for each users
}
  \item{seqfile}{
character, the sequence name you want to analysis
}
  \item{new_project}{
character, the project name
}
}
\details{
Users have to specify a project to upload a job to MG-RAST, just new project supported in this current version of our package.
When the submission process has been successfully completed, MG-RAST ID's ("Accession numbers") will be automatically assigned and the data will be removed from your inbox.
}
\value{
return the MGRAST ID of your sequence and the user name
}
\seealso{
\code{\link{listMgrastInbox}}, \code{\link{uploadMgrast}}
}
\examples{
## login into MGRAST for webkey
# login.info <- loginMgrast("yiluheihei","heihei")

##
#submitMgrastJob(login.info,seq,new_project)
}

