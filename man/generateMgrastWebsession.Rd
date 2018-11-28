\name{generateMgrastWebsession}
\alias{generateMgrastWebsession}
\title{
Generate your websession
}
\description{
Return your websession while login to MG-RAST
}
\usage{
generateMgrastWebsession(cookie)
}
\arguments{
  \item{cookie}{
web cookie, a small piece of data sent from a website and stored in a user's web browser while the user login to MG-RAST
}
}
\details{
A session token is a unique identifier that is generated and sent from a server to a client to identify the current interaction session. The client usually stores and sends the token as an HTTP cookie. This function helps users to get websession from the cookie.
}
\value{
character
}

\seealso{
\code{\link{loginMgrast}}
}
\examples{
##
# generateWebsession(cookie)
}
