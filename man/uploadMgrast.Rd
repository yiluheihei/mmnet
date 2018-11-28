\name{uploadMgrast}
\alias{uploadMgrast}
\title{
Upload your sequence to MGRAST
}
\description{
Sequence can be uploaded to MGRAST for functional annotation.
}
\usage{
uploadMgrast(login.info, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{login.info}{a list of login info generated after user login into MGRAST which consists of webkey of the file, the invlid time of your webkey, cookie, websession and curlhandle. It is unique for each users.
}
  \item{file}{
the sequence file user want to upload, more description in details.
}
}
\details{
MG-RAST supports each user a temporary storage location. This inbox provides temporary storage for data to be submitted for analysis. 
Users can also delete the files in your inbox.

Ensure the file not exsited in the inbox before upload. Uploaded files may be removed from your inbox after 72 hours. Please perform annotation of your files within that time frame. There is a delay between upload completion and appearing in this table due to sequence statistics calculations. This may be on the order of seconds to hours depending on file size. Uses can check your file state with the function \code{\link{listMgrastInbox}}.

Your uploaded sequence data can be in FASTA, FASTQ or SFF format. These are recognized by the file name extension with valid extensions for the appropriate formats .fasta, .fna, .fastq, .fq, and .sff and FASTA and FASTQ files need to be in plain text ASCII. Compressing large files will reduce the upload time and the chances of a failed upload, you can use Zip (.zip) and gzip (.gz) as well as tarred gzipped files (.tgz) but not rar.
We suggest you upload raw data (in FASTQ or SFF format) and let MG-RAST perform the quality control step, see here for details.
}
\seealso{
\code{\link{listMgrastInbox}},\code{\link{submitMgrastJob}}
}
\examples{
## login into MGRAST for webkey
# login.info <- loginMgrast("yiluheihei","heihei")
 
## upload
# uploadMgrast(login.info, file="")

}
