saveMetabolicData <- function(RefDbcache, path = Sys.getenv("HOME")) {
    RefDbcache <- as.list(RefDbcache)
    RefDbcache$user <- Sys.getenv("USERNAME")
    RefDbcache$date <- date()
    RefDbcache$version <- R.Version()
    RefDbcache$date <- date()
    RefDbcache$network <- constructMetabolicNetwork(RefDbcache, path)
    message("saving data to the Specified dir...", domain = NA)
    path <- sprintf("%s/%s", path, ".mmnet")
    if (!file.exists(path)) 
        dir.create(path, showWarnings = FALSE)
    save(RefDbcache, file = file.path(path, "RefDbcache.rda"))
} 
