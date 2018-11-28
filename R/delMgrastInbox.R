delMgrastInbox <- function(login.info, sequence) {
    websession <- login.info$session
    inbox.file <- listMgrastInbox(login.info)
    if (!sequence %in% inbox.file$files){
      warning(sequence," is not in your inbox")
      return(FALSE)
    }
    delete.file <- tryCatch(getForm("http://metagenomics.anl.gov/upload.cgi/user_inbox/?callback=1", 
              auth = websession, websession = "faction", faction = "del", del = "fn", fn = sequence),
              error = function(e) {
                msg <- conditionMessage(e)
                structure(msg, class = "try-error")
              }
    )
    if (inherits(delete.file,"try-error")){
      warning(delete.file)
      return(FALSE)
    }else{
      ret <- delete.file
      str_sub(ret, -3, -1) <- ""
      str_sub(ret, end = 27) <- ""
      ret <- fromJSON(ret)
      if (length(ret$popup_messages) > 0) 
        print(ret$popup_messages)
      return(!(sequence %in% names(ret$locks)))
    }
}
