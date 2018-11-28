uploadMgrast <- function(login.info, file) {
  
  ## set progressbar
  progressDown=function(down, up, pcur, width){
    total=as.numeric(up[1]) # Total size as passed from curlPerform
    cur=as.numeric(up[2])   # Current size as passed from curlPerform
    x=cur/total
    px= round(100 * x)
    ## if(!is.nan(x) &&  px>60) return(pcur) # Just to debug at 60%
    if(!is.nan(x) && px!=pcur){
      x= round(width * x)
      sc=rev(which(total> c(1024^0, 1024^1, 1024^2, 1024^3)))[1]-1
      lb=c('B', 'KB', 'MB', 'GB')[sc+1]
      cat(paste(c("\r  |", rep.int(".", x), rep.int(" ", width - x),
        sprintf("| %g%s of %g%s %3d%%",round(cur/1024^sc, 2), lb, round(total/1024^sc, 2), lb, px)),
        collapse = ""))
      flush.console() # if the outptut is buffered, it will go immediately to console
      return(px)
    }
    return(pcur)
  }
  if (!file.exists(file)) {
      message("Can not find file, please check if file exsits!")
      return(FALSE)
  }
  if (any(grepl(file_path_sans_ext(basename(file)),listMgrastInbox(login.info)$files)))
      stop("sequence file have been exits in your Inbox")
  width <- 50
  pcur  <- 0
  web.feedback <- tryCatch( 
      postForm("http://api.metagenomics.anl.gov/1/inbox/upload", 
        .params = list(upload = fileUpload(file)), 
        .opts = list(httpheader = c(auth = login.info$webkey),noprogress = FALSE,
          progressfunction=function(down,up)
            pcur <<- progressDown(down, up, pcur, width))),
      error = function(e) {
        msg <- conditionMessage(e)
        structure(msg, class = "try-error")
      }
  )
  if (inherits(web.feedback,"try-error")){
    warning(web.feedback)
    return(FALSE)
  }else{
    if (grepl("data\\s+received\\s+successfully", web.feedback)) {
      cat(paste0("\n", file, " successfully upload ", date(), "\n"))
      return(TRUE)
    } else {
      cat("please check your webkey and the file, and ensure uploading file does not exist in your inbox\n")
      return(FALSE)
    }
  }
} 