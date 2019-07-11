library(KEGGREST)
library(igraph)
library(XML)
library(plyr)

# This function preprocess the KO information to delete the redundant metabolites 
preprocessKOMetabolites <- function(ko.info) {
  ko.info <- ko.info[sapply(ko.info, class) != "character"]
  ko.info <- do.call(rbind, ko.info)[c("name", "type", "substrate.name", "product.name")]
  ko.info[, 1] <- gsub("ko:", "", ko.info[, 1])
  ko <- ko.info$name
  multi.ko.info <- gregexpr("(\\s)", ko, perl = TRUE)
  multi.ko.index <- which((lapply(ko, nchar) != 6))
  multi.ko <- strsplit(ko[multi.ko.index], "\\s")
  ko.info[multi.ko.index, 1] <- sapply(multi.ko, function(x) x[1])
  multi.ko <- sapply(multi.ko, function(x) x[-1])
  df <- ko.info[multi.ko.index, -1]
  #multi.rn.metabolites <- rep(ko.info[multi.ko.index, -1], times = listLen(multi.ko))
  multi.rn.metabolites <- df[rep(1:nrow(df), times = Biobase::listLen(multi.ko)),]
  multi.ko.info <- cbind(unlist(multi.ko), multi.rn.metabolites)
  names(multi.ko.info) <- names(ko.info)
  ko.info <- rbind(ko.info, multi.ko.info)
  message("Processing of reversible reaction  ...", domain = NA)
  reverse.index <- which(ko.info$type == "reversible")
  ko.info$substrate.name[reverse.index] <- paste(ko.info$substrate.name[reverse.index], 
    ko.info$product.name[reverse.index], sep = ",")
  ko.info$product.name[reverse.index] <- ko.info$substrate.name[reverse.index]
  message("Merging the substrates and products of the same enzyme  ...", domain = NA)
  delete.no <- which(duplicated(ko.info$name))
  dup <- unique(ko.info$name[-(Biobase::isUnique(ko.info$name))])
  dup.select <- match(dup, ko.info$name)
  dup.no <- sapply(dup, function(x) which(ko.info$name == x))
  dup.no <- sapply(dup.no, function(x) x[-1])
  ko.info[dup.select, 3] <- mapply(function(x, y) paste(c(ko.info[3][x, ], ko.info[3][y, 
    ]), collapse = ","), dup.select, dup.no)
  ko.info[dup.select, 4] <- mapply(function(x, y) paste(c(ko.info[4][x, ], ko.info[4][y, 
    ]), collapse = ","), dup.select, dup.no)
  ko.info <- ko.info[-delete.no, ]
  message("Deleting the same metabolites for each enzyme  ...", domain = NA)
  dlt.dup.metabolites <- function(x) {
    str.metabio <- strsplit(x, ",")
    return(unique(unlist(str.metabio)))
  }
  ko.info$substrate.name <- sapply(ko.info[, 3], dlt.dup.metabolites)
  ko.info$product.name <- sapply(ko.info[, 4], dlt.dup.metabolites)
  names(ko.info) <- c("ko", "type", "substrate", "product")
  delete.ko.no <- which(lapply(ko.info$ko, nchar) != 6)
  if (length(delete.ko.no)) {
    return(ko.info[-delete.ko.no, -2]) 
  } else { 
    return(ko.info[, -2])
  }
} 

# construct global reference network
constructMetabolicNetwork <- function(RefDbcache) {
  product <- RefDbcache$product
  substrate <- RefDbcache$substrate
  message("construct the edge matrix of ko according to metabolites")
  edge.matrix <- matrix(0, length(product), length(product))
  edge.matrix <- laply(product, function(y) laply(substrate, function(x) length(intersect(y, 
    x))), .progress = "text")
  edge.matrix[edge.matrix > 1] <- 1
  diag(edge.matrix) <- 0
  rownames(edge.matrix) <- colnames(edge.matrix) <- RefDbcache$ko
  message("constrcuting the ref Network ")
  g <- graph.adjacency(edge.matrix, mode = "directed")
  g <- set.graph.attribute(g, name = "name", value = "refNet")
  return(g)
} 

# extract kegg pathway info
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
  } else {
    ko.info <- "No reaction in this pathway"
  }
  
  ko.info
} 


ko.list <- keggList("pathway/ko")
deleteGlobal.info <- regexpr("ko01|ko00", names(ko.list), perl = TRUE)
ko.list <- ko.list[deleteGlobal.info > 0]
ko.path <- regexpr("ko\\d+", names(ko.list), perl = TRUE)
ko.pathway <- substring(names(ko.list), ko.path, ko.path + attr(ko.path, "match.length") - 
    1)
RefDbcache <- lapply(ko.pathway, getKOPathwayInfo)
names(RefDbcache) <- ko.pathway
RefDbcache <- as.list(preprocessKOMetabolites(RefDbcache))
RefDbcache <- as.list(RefDbcache)
RefDbcache$user <- Sys.getenv("USERNAME")
RefDbcache$date <- date()
RefDbcache$version <- R.Version()
RefDbcache$date <- date()
RefDbcache$network <- constructMetabolicNetwork(RefDbcache)
message("saving data to the Specified dir...", domain = NA)
# path <- sprintf("%s/%s", path, ".mmnet")
save(RefDbcache, file = "../R/sysdata.rda")
