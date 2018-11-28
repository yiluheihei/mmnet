loadMetabolicData <- function(path = Sys.getenv("HOME")) {
    path <- sprintf("%s/%s", path, ".mmnet")
    if (!exists(path)) {
        message(" Attempt to load default data.\nRecommand run updateKEGGPathway() to get the latest version of data."
                , domain = NA)
        data(RefDbcache)
    } else load(file = sprintf("%s/%s", path, "RefDbcache.rda"))
} 
