\name{generateMgrastWebkey}
\alias{generateMgrastWebkey}
\title{
Generate your webkey on MGRAST
}
\description{
This function help users to generate a webkey which is essential for data uploding to MGRAST.
}
\usage{
generateMgrastWebkey(login.info)
}
\arguments{
  \item{login.info}{
list, reponse information while user login to MGRAST, contains webkey, invalid time of webkey, cookie, websseion and curlhandle, see \code{\link{loginMgrast}} for more details}
}
\value{
\item{webkey}{character, it will be valid for a limited time only (the next value)}
\item{invalid.until}{character, with corresponding to the webkey, just the limit time of your webkey}
}
\seealso{
\code{\link{loginMgrast}}
}
\examples{
# not run
# generateMgrastWebkey(login.info$curlhandle)
}
