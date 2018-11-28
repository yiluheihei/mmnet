\name{listMgrastInbox}
\alias{listMgrastInbox}
\title{
List the contents of the user inbox 
}
\description{
List the file and file information user uploaded to MGRAST.
}
\usage{
listMgrastInbox(login.info)
}
\arguments{
  \item{login.info}{
list, reponse information while user login to MGRAST, contains webkey, invalid time of webkey, cookie, websseion and curlhandle, see \code{\link{loginMgrast}} for more details}
}
\details{
There is a delay between upload completion and before job submission because sequence statistics is begin to calculate once file uploaded . This may be on the order of seconds to hours depending on file size. This files cannot be submitted to create job and annotation until analysis is complete. File name, file (calculation completed) statistics information and under analysis files which appear in the locks will be list with this function.
}
\value{
a list, contains file details user uploaded
}


\seealso{
\code{\link{submitMgrastJob}}
}
\examples{
## login into MGRAST for webkey
# login.info <- loginMgrast("yiluheihei","heihei")

## not run need user login to MGRAT
#listMgrastInbox(login.info)
}
