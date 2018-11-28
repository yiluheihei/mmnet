loginMgrast <- function(user, userpwd) {
    curlGlobalInit()
    if (class(user) != "character" || class(userpwd) != "character") 
        stop("both user and user password must be character")
    http_header <- c(Accept = "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8", 
        `Accept-Encoding` = "gzip,deflate,sdch", `Accept-Language` = "zh-CN,zh;q=0.8", 
        Connection = "keep-alive", `Content-Type` = "multipart/form-data; boundary=----WebKitFormBoundaryrcMPojVmeWEDcGX8", 
        `User-Agent` = "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/30.0.1599.101 Safari/537.36")
    ch <- getCurlHandle()
    data_login <- c(page = "Login", action = "perform_login", login = user, password = userpwd)
    # tryCatch(login <- postForm(uri = "http://metagenomics.anl.gov/", .params = data_login, 
    #                 .opts = list(cookiejar = "cookie"), curl = ch, httpheader = http_header, verbose = TRUE),
    #          error = function(err) {stop(simpleError("Your Internet not connected or MGRAST host can not be connetecd, please try later"))}
    # )
    login <- tryCatch(postForm(uri = "http://metagenomics.anl.gov/", .params = data_login, 
                   .opts = list(cookiejar = "cookie"), curl = ch, httpheader = http_header, verbose = TRUE),
             error = function(e) {
                msg <- conditionMessage(e)
                structure(msg, class = "try-error")
             }
     )
    if (inherits(login,"try-error")){
      warning(login)
      return(NULL)
    }else{
        cookie <- getCurlInfo(ch)[["cookielist"]]
        if (grepl("page=AccountManagement", login, perl = TRUE)) {
            cat(paste0(user, " login successed\n"))
            ret <- generateMgrastWebkey(ch)
            ret$cookie <- cookie
            ret$session <- strsplit(ret$cookie, "\t")[[1]][[7]]
            ret$curlhandle <- ch
            return(ret)
        } else {
            cat("Login failed, user name or password error\n")
            return(NULL)
        }
    }
} 
