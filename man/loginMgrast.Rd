\name{loginMgrast}
\alias{loginMgrast}
\title{
Login to the MGRAST
}
\description{
Login to MGRAST, obtain your Websession and webkey for data uploading, as an authenticate of subsequent operation about MGRAST (e.g. data upload, job submission). 
}
\usage{
loginMgrast(user, userpwd)
}
\arguments{
  \item{user}{
character, username you have registered on MGRAST

}
  \item{userpwd}{
character, password correspondance to your username
}

}
\details{
Sequence functional annotation is the most necessary and important step for metagenomic analysis. MGRAST (metagenomics RAST) server is a free and public utomated analysis platform for metagenomes providing quantitative insights into microbial populations based on sequence data, also very fast, stable and extensible. Thus, MGRAST was took for sequence functional annotation in our package \emph{mmnet}.

Webkey, session and curlhandle would be generated when users login to MGRAST with this function, which is the key for data uploading and sequence annotation.
}
\value{
\item{webkey}{character, e.g. "8Dvg2d5DCp7KsWKBPzY2GS4i", the webkey of your account which is essential for data uploding to MGRAST. However, it will be valid for a limited time only (the next value). Users will generate a new webkey with fucntion \code{\link{generateMgrastWebkey}} for convenient.}
\item{invalid.until}{character, with corresponding to the webkey, just the limit time of your webkey}
\item{cookie}{data saved on user's computer while login to MGRAST }
\item{session}{character, a userid that MGRAST server assigned}
\item{curlhandle}{class of CURLHandle, the curl handle while users to MGRAST which is saved for user's sebsequent analysis with MGRAST.}
}
\seealso{
\code{\link{uploadMgrast}}
}
\examples{
# need a username and password on MGRAST
# login.info <- loginMgrast("mmnet","mmnet")
}

