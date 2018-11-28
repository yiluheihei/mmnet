\name{delMgrastInbox}
\alias{delMgrastInbox}
\title{
Delete the files in you MG-RAST inbox
}
\description{
Once job is submiited to MG-RAST, user can delete the redundant sequence.
}
\usage{
delMgrastInbox(login.info, sequence)
}
\arguments{
  \item{login.info}{
web information while login to MGRAST , more details see \code{\link{loginMgrast}}  
}
  \item{sequence}{
character, represents the file your will delete 
}
}
\details{
Users can delete files from the inbox on MG-RAST, Filenames in your Inbox are undergoing analysis and cannot be moved or submitted to a different process until analysis is complete.
}

