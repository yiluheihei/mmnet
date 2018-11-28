constructMetabolicNetwork <- function(path = Sys.getenv("HOME")) {
    loadMetabolicData(path)
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
