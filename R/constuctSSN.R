constructSSN <- function(abundance) {
    if (!exists("RefDbcache"))
         loadMetabolicData(path = Sys.getenv("HOME"))
    RefDbcache <- get("RefDbcache", envir = parent.frame())
    if (!is.igraph(RefDbcache$network)) 
        stop("not a igraph object")
    if (!inherits(abundance, "biom")) 
        stop("abundance must be BIOM format")
    refnode <- V(RefDbcache$network)$name
    if(abundance$matrix_type == "sparse"){
      abundance.data = as.matrix(biom_data(abundance))
    }else{
      abundance.data = as.matrix(biom_data(abundance)) 
    }
    single_g <- function(biom.data, abund.data){
      subnodes <- intersect(unlist(observation_metadata(biom.data))[abund.data != 0], refnode)
      if (!length(subnodes)) 
        stop("names of abundance should be KO number or there is no KO intersection between this sample and reference data")
      g <- induced.subgraph(RefDbcache$network, subnodes)
      match.index <- match(V(g)$name, unlist(observation_metadata(biom.data)))
      subabund <- abund.data[match.index]
      g <- set.graph.attribute(g, "name", "SSN")
      g <- set.vertex.attribute(g, "abundance", index = V(g), value = subabund)
      g <- delete.vertices(g, names(which(igraph::degree(g, mode = "all") == 0)))
      return(g)
    }
    if (ncol(abundance) == 1 ){
      g <- single_g(abundance,abundance.data)
    }else{
      g <- apply(abundance.data, 2, function(x)single_g(abundance,x))
    }
    return(g)
} 
