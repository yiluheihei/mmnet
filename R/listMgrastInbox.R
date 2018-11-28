listMgrastInbox <- function(login.info) {
    websession <- login.info$session
    file.inform <- tryCatch(
      getForm("http://metagenomics.anl.gov/upload.cgi/user_inbox/?callback=1", 
              auth = websession),
      error = function(e) {
        msg <- conditionMessage(e)
        structure(msg, class = "try-error")
      }
    )
    if (inherits(file.inform, "try-error")){
      warning(file.inform)
      return(FALSE)
    }else{
      if (length(grep("unauthorized\\s+request", file.inform)) != 0) {
        stop("the authorize may incorrect")
      }
      str_sub(file.inform, -3, -1) <- ""
      str_sub(file.inform, end = 27) <- ""
      return(fromJSON(file.inform))
    }
} 
