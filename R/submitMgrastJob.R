submitMgrastJob <- function(login.info, seqfile, new_project) {
    file <- listMgrastInbox(login.info)
    if (!match(seqfile, file$files, nomatch = 0)){
      if (!match(file_path_sans_ext(seqfile), file$files, nomatch = 0)){
        stop(paste0(seqfile, " is not in your Mgrast Inbox"))
        return(FALSE)
      }
    }
    ## unpack the sequence file
    websession <- login.info$session
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
    if (grepl("(\\.gz$)|(\\.zip$)|(\\.tgz$)", seqfile))
    {
        message("Unpack the compressed File...",domain = NULL)
        while (seqfile %in% file$files){
          Sys.sleep(10)
          unpackfile <- tryCatch(getForm("http://metagenomics.anl.gov/upload.cgi/user_inbox/?callback=1", 
            auth = websession, websession = "faction", faction = "unpack", unpack = "fn", fn=seqfile),
            error = function(e) {
              msg <- conditionMessage(e)
              structure(msg, class = "try-error")
            }
          )
          if (inherits(unpackfile, "try-error")){
            warning(unpackfile)
            return(FALSE)
          }else{
            file = listMgrastInbox(login.info)
          }
        }
        seqfile <- file_path_sans_ext(seqfile)
    }

    message("Under statistics caculation...", domain = NULL)
    while (seqfile %in% names(file$locks)) {
      # stop('The fn is under statistics calculation')
      refresh <- tryCatch(getForm("http://metagenomics.anl.gov/upload.cgi/user_inbox/?callback=1", 
                                  auth = websession, websession = "faction"),
                          error = function(e) {
                            msg <- conditionMessage(e)
                            structure(msg, class = "try-error")
                          }
      )
      if (inherits(refresh, "try-error")){
        warning(refresh)
        return(FALSE)
      }else{                   
        Sys.sleep(10)
        file <- listMgrastInbox(login.info)
      }
    }
    message("Statistics caculation is complete...", domain = NULL)
    check.dup <- tryCatch(getForm("http://metagenomics.anl.gov/metagenomics.cgi?page=Upload",auth=login.info$session,
        action="check_for_duplicates",curl = login.info$curlhandle, seqfiles=seqfile), 
        error = function(e) {
          msg <- conditionMessage(e)
          structure(msg, class = "try-error")
        }
    )
    if (inherits(check.dup, "try-error")){
      warning(check.dup)
      return(FALSE)
    }else{
      if (grepl("already exist in MG-RAST",check.dup)){
        message(paste(seqfile, "already exist in MG-RAST","\n", "existing ID:", str_extract(check.dup, "\\d+\\.\\d+")))
        metagenome.id <- str_extract(check.dup, "\\d+\\.\\d+")
        names(metagenome.id) <- "metagenome.id"
        return(metagenome.id)
      }
      dereplication <- "dereplication"
      screening <- "h_sapiens_asm"
      filter_ln <- "filter_ln"
      deviation <- 2
      filter_ambig <- "filter_amibg"
      max_ambig <- 5
      # seqfile <- file_path_sans_ext(seqfile)
      if (missing(new_project)) 
        new_project <- as.character(Sys.time())
      para <- c(create_job = 1, table_perpage_0 = 10, seqfiles = seqfile, dereplication = dereplication, 
                screening = screening, filter_ln = filter_ln, deviation = deviation, filter_ambig = filter_ambig, 
                max_ambig = max_ambig, priorityOption = "never")
      para <- c(para, new_project = new_project)
      message("Create the MGRAST project...", domain= NA )
      response <- tryCatch(getForm("http://metagenomics.anl.gov/metagenomics.cgi?page=Upload", 
                                   .params = para, curl = login.info$curlhandle),
                                    error = function(e) {
                                      msg <- conditionMessage(e)
                                      structure(msg, class = "try-error")
                                    }
      )
      if (inherits(response, "try-error")){
        warning(response)
        return(FALSE)
      }else{
        response <- htmlParse(response)
        user <- getNodeSet(response, "//div[@id='user']")
        account <- sapply(user, xmlValue)
        account <- gsub("\\n+", "", account)
        metagenome.id <- getNodeSet(response, "//div[@class='modal-body']")
        metagenome <- sapply(metagenome.id, function(x) sapply(x["p"], xmlValue))
        metagenome <- unlist(try(metagenome[which(grepl("success", metagenome))], silent = T))
        names(metagenome) <- NULL
        if (length(metagenome)) {
          print(metagenome[1])
        } else {
          stop("job submmission failed")
        }
        metagenome.id <- str_extract(metagenome[2], "\\d+\\.\\d+")
        names(metagenome.id) <- "metagenome.id"
        return(metagenome.id)
        message("Job submiited", domain= NA )
      }
    }
    
} 
