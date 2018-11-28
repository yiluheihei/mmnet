generateMgrastWebsession <- function(cookie) {
    cookie <- as.matrix(read.delim(cookie, header = F, stringsAsFactors = F))
    key.index <- match("WebSession", cookie)
    if (key.index) {
        key <- cookie[key.index + nrow(cookie)]
    } else {
        stop("cookie may incorrect")
    }
    return(key)
} 
