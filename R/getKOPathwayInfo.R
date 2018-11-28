getKOPathwayInfo <- function(koPathway) {
    pathway.kgml <- keggGet(koPathway, "kgml")
    pathway.kgml <- xmlRoot(xmlTreeParse(pathway.kgml))
    attr <- xmlSApply(pathway.kgml, xmlAttrs)
    reaction <- as.matrix(attr[names(attr) == "reaction"])
    if (length(reaction)) {
        reaction <- do.call(rbind, reaction)
        entry <- attr[names(attr) == "entry"]
        entry <- na.omit(t(sapply(entry, function(x) x[c("id", "name", "reaction")])))
        relation <- attr[names(attr) == "relation"]
        metabolites <- lapply(pathway.kgml[names(pathway.kgml) == "reaction"], xmlChildren)
        metabolites <- lapply(metabolites, function(x) lapply(x, xmlAttrs, "name"))
        reaction.info <- lapply(lapply(lapply(metabolites, unlist), function(x) tapply(x, 
            names(x), c)), function(x) return(cbind(paste(x[["substrate.name"]], 
            collapse = ","), paste(x[["product.name"]], collapse = ","))))
        reaction.info <- do.call(rbind, reaction.info)
        reaction <- cbind(reaction, reaction.info)
        colnames(reaction) <- c("id", "rn", "type", "substrate.name", "product.name")
        entry <- entry[is.element(entry[, 1], reaction[, 1]), ]
        colnames(entry) <- c("id", "name", "rn")
        reaction <- apply(reaction, 2, function(x) x[match(entry[, 1], reaction[, 
            1])])
        ko.info <- as.data.frame(cbind(entry, reaction), row.names = F, stringsAsFactors = FALSE)
        message(paste(koPathway, "processing is complete"))
        ko.info <- ko.info[unique(names(ko.info))]
    } else ko.info <- "No reaction in this pathway"
    return(ko.info)
    if (FALSE) {
        a <- getNodeSet(xmlParse(ko00010), "//reaction")
        xmlSApply(ko, xmlValue)
        getChildrenStrings(getNodeSet(xmlParse(ko00010), "//reaction"))
        sapply(getNodeSet(xmlParse(ko00010), "//reaction"), xmlAttrs)
        sapply(getNodeSet(xmlParse(ko00010), "//reaction"), function(x) xmlAttrs(getSibling(x)))
        ko.list <- keggList("pathway/ko")
        deleteGlobal.info <- regexpr("ko01|ko00", names(ko.list), perl = TRUE)
        ko.list <- ko.list[deleteGlobal.info > 0]
        ko.info <- regexpr("ko\\d+", names(ko.list), perl = TRUE)
        ko.pathway <- substring(names(ko.list), ko.info, ko.info + attr(ko.info, 
            "match.length") - 1)
    }
} 
