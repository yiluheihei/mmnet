checkMgrastMetagenome <- function(login.info, metagenome.id, public = TRUE){
  if (public)
    url <- paste0("http://api.metagenomics.anl.gov/1/metagenome/", paste0("mgm",
                metagenome.id), "?verbosity=metadata")
  else{
    if (missing(login.info)){
      warning("Private project checking need users login to MGRAST")
      return(NULL)
    }
    url <- paste0("http://api.metagenomics.anl.gov/1/metagenome/", paste0("mgm",
                metagenome.id), "?verbosity=metadata", "&auth=", login.info$webkey)
  }
  url.response <- tryCatch(getURL(url),
                           error = function(e) {
                             msg <- conditionMessage(e)
                             structure(msg, class = "try-error")
                           }
  )
  if (inherits(url.response, "try-error")){
    warning(url.response)
    return(FALSE)
  }else{
    if (grepl("ERROR", url.response)){
      message(url.response)
      return(metagenome.id = FALSE)
    }else{
      return(metagenome.id = TRUE)
    }
  }
}
