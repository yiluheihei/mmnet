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
    multi.rn.metabolites <- df[rep(1:nrow(df), times = listLen(multi.ko)),]
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
    dup <- unique(ko.info$name[-(isUnique(ko.info$name))])
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
    if (delete.ko.no) 
        return(ko.info[-delete.ko.no, -2]) else return(ko.info[, -2])
} 