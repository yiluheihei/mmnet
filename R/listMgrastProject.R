listMgrastProject <- function(login.info) {
    response <- tryCatch(
      getForm("http://metagenomics.anl.gov/metagenomics.cgi", page = "Upload", 
        curl = login.info$curlhandle),
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
      user.project <- xmlDoc(getNodeSet(response, "//div[@id='sel_project_div']")[[1]])
      job.project <- xpathApply(user.project, "//option")[-1]
      if (!length(job.project)) {
        stop("You have never submmitted a job in Mgrast\n")
      }
      job.project <- sapply(job.project, toString.XMLNode)
      greg <- sapply(job.project, function(x) gregexpr("value=\\\"(\\d+)\\\">(.*)<", 
                                                       x, perl = T))
      job <- mapply(function(x, y) mapply(function(start, length) substring(x, start, 
                                                                            start + length - 1), attr(y, "capture.start"), attr(y, "capture.length")), 
                    job.project, greg)
      colnames(job) <- job[2, ]
      project.url <- "http://api.metagenomics.anl.gov/1/project/mgp"
      project.url <- sapply(job[1, ], function(x) paste0(project.url, x, "?verbosity=full&auth=", 
                                                         login.info$webkey))
      project <- lapply(lapply(project.url, getURL), fromJSON)
      project <- lapply(project, function(x) x[c("name", "metagenomes", "id")])
      return(project)
    }    
} 
