updateKEGGPathway <- function(path = Sys.getenv("HOME")) {
    ko.list <- keggList("pathway/ko")
    deleteGlobal.info <- regexpr("ko01|ko00", names(ko.list), perl = TRUE)
    ko.list <- ko.list[deleteGlobal.info > 0]
    ko.path <- regexpr("ko\\d+", names(ko.list), perl = TRUE)
    ko.pathway <- substring(names(ko.list), ko.path, ko.path + attr(ko.path, "match.length") - 
        1)
    RefDbcache <- lapply(ko.pathway, getKOPathwayInfo)
    names(RefDbcache) <- ko.pathway
    RefDbcache <- as.list(preprocessKOMetabolites(RefDbcache))
    saveMetabolicData(RefDbcache, path)
} 
