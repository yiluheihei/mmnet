topologicalAnalyzeNet <- function(g, Scatterplot = TRUE, .properties = c("betweennessCentrality", 
                        "degree", "clusteringCoefficient", "pageRank"), 
                        mode = c("all", "in", "out"), ...) {
  
  ## plot scatter matrix
  plotScattermat <- function(df){
      df[df <= 0] <- NA
      dat <- melt(df, names(df)[1])
      #dat [dat == 0] <- NA
      nan.index <- sapply(dat[c(1,3)], is.nan)
      infinite.index <- sapply(dat[c(1,3)], is.infinite)
      na.index <- is.na(dat[c(1,3)])
      invalid.index <- !nan.index & !infinite.index & !na.index
      invalid.index <- invalid.index[,1] & invalid.index[,2]
      dat <- dat[invalid.index,] 
      print(ggplot(dat, aes_string(x ="value", y = names(df)[1])) 
        + geom_point(na.rm = TRUE) + stat_smooth(method = "lm", colour = "red", na.rm = TRUE) + facet_wrap(~variable, 
          ncol = 2, scales = "free") + scale_x_log10() +  scale_y_log10())
  }


  if (!is.igraph(g)) {
    stop("Not a igraph object")
  }
  .properties <- c(.properties, c(...))
  if (length(.properties) == 0) 
    stop("No inputs passed to properties")
    proper.options <- c("betweennessCentrality", "degree", "clusteringCoefficient", 
                        "pageRank")
    property <- match(proper.options, unlist(.properties), nomatch = 0)
    if (max(property) == 0) 
      stop("input is invalid", domain = NULL)
    mode <- match.arg(mode, c("all", "in", "out"))
  vertex.attr <- list.vertex.attributes(g)

  if (! "abundance" %in%  vertex.attr){
    vertex.attribute <- list.vertex.attributes(g)
    if (match("p.value", vertex.attribute))
      diff.attribute   <- setdiff(vertex.attribute, c("name", "p.value"))
    else
      diff.attribute <- setdiff(vertex.attribute, "name")
    diff.abundance <- get.vertex.attribute(g, name = diff.attribute, index = V(g))
    topo.params <- data.frame(matrix(0, length(diff.abundance), length(.properties) +1))
    names(topo.params) <- c(diff.attribute, unlist(.properties))
    topo.params[diff.attribute] <- diff.abundance
  }else{
      abundance <- get.vertex.attribute(g, name = "abundance", index = V(g))
      topo.params <- data.frame(matrix(0, length(abundance), length(.properties) +1))
      names(topo.params) <- c("abundance", unlist(.properties))
      topo.params[1] <- abundance
  }
  ftopo.value <- function(x, type) {
  switch(type, betweennessCentrality = igraph::betweenness(x), 
          degree = igraph::degree(x, mode = mode),
          clusteringCoefficient = igraph::transitivity(x, type ="local", 
          vids = V(x), isolates=c("zero")),
          pageRank  = igraph::page.rank(x)$vector)
  }            
  for (para in .properties){
    if (match(para, proper.options, nomatch = 0)){
      topo.params[para] <- ftopo.value(g, para)
      g <- set.vertex.attribute(g, name = para, index = V(g), value = topo.params[para])
    }
  } 
  if (Scatterplot) 
    plotScattermat(topo.params)
  return(g)

}

    