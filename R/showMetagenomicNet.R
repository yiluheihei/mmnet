showMetagenomicNet <- function(net, mode = c("ref", "ssn", "compared"), method = c("OR", "rank", 
    "JSD"), cutoff, ...) {
    ## vertex color for compared network
    VertexColor <- function(value, method, cutoff) {
        vertex.color <- rep("grey", length(value))
        if (method == "OR") {
            value <- (log2(value))
            enrich.no <- which(value > log2(cutoff[[2]]))
            deplete.no <- which(value < log2(cutoff[[1]]))
        } else {
            deplete.no <- which(value < quantile(value, probs = cutoff[[1]]))
            enrich.no <- which(value > quantile(value, probs = cutoff[[2]]))
        }
        vertex.color[enrich.no] <- "red"
        vertex.color[deplete.no] <- "green"
        return(vertex.color)
    }
    
    net <- delete.vertices(net, names(which(igraph::degree(net) == 0)))
    method <- match.arg(method, c("OR", "rank", "JSD"))
    if (!is.igraph(net)) 
        stop("not a igraph object")
    mode <- match.arg(mode, c("ref", "ssn", "compared"))
    if (mode == "ssn") {
        tmp <- get.vertex.attribute(net, name = "abundance", index = V(net)) + 1
        if (is.null(tmp)) 
            stop("error mode, please check your graph mode")
        vertex.size2 <- c(1:length(tmp))
        vertex.size2[which(tmp == 1)] <- 2
        vertex.size2[which(tmp != 1)] <- 2 + log(tmp[which(tmp != 1)])
        plot.igraph(net, vertex.size = vertex.size2, ...)
    }
    if (mode == "ref") {
        plot.igraph(net, ...)
    }
    if (mode == "compared") {
        if (method == "OR") {
            OR <- get.vertex.attribute(net, name = "OR", index = V(net))
            if (is.null(OR)) 
                stop("error mode, please check your graph mode")
            vertex.color <- VertexColor(OR, method = "OR", cutoff = c(0.5, 2))
            plot.igraph(net, vertex.color = vertex.color, ...)
        }
        if (method == "rank") {
            diffabund <- get.vertex.attribute(net, name = "diffabund", index = V(net))
            if (is.null(diffabund)) 
                stop("error mode, please check your graph mode")
            vertex.color <- VertexColor(diffabund, method = "rank", cutoff = c(0.1, 
                0.9))
            plot.igraph(net, vertex.size = 3, vertex.color = vertex.color, ...)
        }
        if (method == "JSD") {
            JSD <- get.vertex.attribute(net, name = "JSD", index = V(net))
            if (is.null(JSD)) 
                stop("error mode, please check your graph mode")
            vertex.color <- VertexColor(JSD, method = "JSD", cutoff = c(0.1, 0.9))
            if (is.null(JSD)) 
                stop("error mode, please check your graph mode")
            plot.igraph(net, vertex.color = vertex.color, ...)
        }
    }
    
} 
