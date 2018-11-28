## ----include=FALSE--------------------------------------------
library(knitr)
options(width=64,digits=2)
opts_chunk$set(size="small")
opts_chunk$set(tidy=TRUE,tidy.opts=list(width.cutoff=50,keep.blank.line=TRUE))
opts_knit$set(eval.after='fig.cap')
# for a package vignette, you do want to echo.
# opts_chunk$set(echo=FALSE,warning=FALSE,message=FALSE)
opts_chunk$set(warning=FALSE,message=FALSE)
opts_chunk$set(cache=TRUE,cache.path="cache/mmnet")

## ----install-pkg, eval=FALSE----------------------------------
#  ## install release version of mmnet
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("mmnet")
#  
#  ##install the latest development version
#  useDevel()
#  biocLite("mmnet")

## ----load-pkg,eval=TRUE, include=FALSE------------------------
library(mmnet)

## ----load-pkg2,eval=FALSE-------------------------------------
#  library(mmnet)

## ----sample-download, eval=FALSE------------------------------
#  download.file("ftp://ftp.metagenomics.anl.gov/projects/10/4440616.3/raw/507.fna.gz",
#                destfile = "obesesample.fna.gz")
#  download.file("ftp://ftp.metagenomics.anl.gov/projects/10/4440823.3/raw/687.fna.gz",
#                destfile = "leansample.fna.gz")

## ----login-mgrast, eval = FALSE, echo=TRUE--------------------
#  ## login on MG-RAST
#   login.info <- loginMgrast(user="mmnet",userpwd="mmnet")

## ----upload-mgrast, eval=FALSE,echo=TRUE----------------------
#  ## select the sample data to upload
#  seq.file <- c("obesesample.fna.gz", "leansample.fna.gz")
#  ## upload sequence data to MG-RAST
#  metagenome.id <- lapply(seq.file, uploadMgrast, login.info = login.info)

## ----submit-mgrast, eval=FALSE, echo=TRUE---------------------
#  ## submit MGRAST job
#    metagenome.id <- lapply(seq.file, submitMgrastJob, login.info = login.info)
#    show(metagenome.id)

## ----check-metagenome, eval=FALSE,echo=TRUE-------------------
#  ## check MGRAST project status
#    metagenome.status <- lapply(metagenome.id, checkMgrastMetagenome, login.info = login.info)
#  ## apparently, status of completed annotation metagenome is TRUE
#    show(metagenome.status)

## ----getannotation,eval=FALSE, echo=TRUE----------------------
#  ## private data
#    private.annotation <- lapply(metagenome.id, getMgrastAnnotation, login.info=login.info)
#  ## public annotation data, does not require login.info
#    public.annotation <- lapply(metagenome.id, getMgrastAnnotation)

## ----data-load------------------------------------------------
data(anno)
summary(anno)

## ----MG-RAST-mmnet, eval=FALSE--------------------------------
#  ## first login on MG-RAST
#  login.info <- loginMgrast("mmnet", "mmnet")
#  ## prepare the metagenomic sequence for annotation
#  seq <- "obesesample.fna.gz"
#  ## mgrast annotation
#  uploadMgrast(login.info, seq)
#  metagenome.id2 <- submitMgrastJob(login.info, seqfile = basename(seq))
#  while (TRUE) {
#    status <- checkMgrastMetagenome(metagenome.id = metagenome.id2)
#    if (status)
#      break
#    Sys.sleep(5)
#    cat("In annotation, please waiting...")
#    flush.console()
#  }
#  ## if annotation profile is public,take login.info as NULL
#  anno2 <- getMgrastAnnotation(metagenome.id2, login.info=login.info)

## ----estimate, echo=TRUE--------------------------------------
 mmnet.abund <- estimateAbundance(anno)
 show(mmnet.abund)

## ----MG-RAST-BIOM, fig.align='center',fig.cap='enzymatic genes abundance comparison between MG-RAST and mmnet package',dev='pdf',fig.show='hold',out.width='.7\\linewidth', out.height='.7\\linewidth'----
## download BIOM functional abundance profile of the two sample metagenome from MG-RAST
if(require(RCurl)){
  function.api <- "http://api.metagenomics.anl.gov/1/matrix/function"
  mgrast.abund <- read_biom(getForm(function.api,
                             .params=list(id="mgm4440616.3",id="mgm4440823.3",
                                          source="KO",result_type="abundance")))
}
## obtain the intersect ko abundance of MG-RAST and esimatiAbundance
intersect.ko <- intersect(rownames(mgrast.abund),rownames(mmnet.abund))
## compare the two by taking one metagenome
mgrast.abund1 <- biom_data(mgrast.abund)[,1][intersect.ko]
mmnet.abund1 <- biom_data(mmnet.abund)[,1][intersect.ko]
if(require(ggplot2)){
  p <- qplot(mgrast.abund1,mmnet.abund1) + 
            geom_abline(slope=1, intercept=0,color="blue") +
            ylim(0, 400) + xlim(0, 400)
  print(p)
}

## ----IMGSample, eval=FALSE------------------------------------
#  ## Load the IMG/M sample data
#  IMGFile <- system.file("extdata/IMGSample.tab",package="mmnet")
#  IMGSample <- read.delim2(IMGFile,quote = "")
#  ## Create BIOM file for network construction
#  abundance <- IMGSample$Gene.Count
#  abundance <- abundance/sum(abundance)
#  abundance <- as.data.frame(abundance)
#  KO <- IMGSample$KO.ID
#  KO <- as.data.frame(gsub("KO:","",KO))
#  biom.data <- make_biom(abundance, observation_metadata = KO)
#  biom.data$type <- "enzymatic genes abundance"

## ----IMGSample2, eval=FALSE-----------------------------------
#  ## Construct and analyze SSN
#  ssn <- constructSSN(biom.data)
#  topologicalAnalyzeNet(ssn)

## ----initial-data,echo=TRUE-----------------------------------
loadMetabolicData()
summary(RefDbcache)

## ----construct-network,echo=TRUE------------------------------
  ssn <- constructSSN(mmnet.abund)
  g <- ssn[[1]]
  summary(g)
  abund <- get.vertex.attribute(g,"abundance",index=V(g))
  summary(abund)

## ----topologicalAnalyzeNet,echo=TRUE,fig.align='center',fig.cap='Topological metabolic network analysis, linking topological features and enzymatic gene abundances',dev='pdf',fig.show='hold',out.width='.7\\linewidth', out.height='.7\\linewidth'----
topo.net <- topologicalAnalyzeNet(g)
## network with topological features as attributes
topo.net

## ----differential,echo=TRUE,fig.align='center',fig.cap='Differential metabolic network analysis, enzymatic genes that are associated with specific state appear as colored nodes (red=enriched, green=depleted)',dev='pdf',fig.show='hold',out.width='.7\\linewidth', out.height='.7\\linewidth'----
state <- c("obese", "lean")
differential.net <- differentialAnalyzeNet(ssn, sample.state= state, method="OR", cutoff = c(0.5, 2))
summary(differential.net)

## ----showMetagenomicNet,echo=TRUE,fig.align='center',fig.cap='Visualization of State Specific Network  using \textit{showMetagenomicNet} with node size  proportional to log(abundance)',dev='pdf',fig.show='hold', out.width='.7\\linewidth', out.height='.7\\linewidth'----
## the reference network
# showMetagenomicNet(RefDbcache$network,mode="ref")
#the state specific metabolic network
showMetagenomicNet(g, mode="ssn", vertex.label = NA, edge.width = 0.3, edge.arrow.size = 0.1, edge.arrow.width = 0.1, layout = layout.fruchterman.reingold)

## ----RCytoscape, eval = FALSE, echo = TRUE--------------------
#  if (require(RCytoscape)){
#    refnet <- ssn[[1]]
#    net <- igraph.to.graphNEL(refnet)
#    ## initialize the edge attibute
#    #edge.attr=list.edge.attributes(refnet)
#    #edge.attr.class = sapply(edge.attr, class)
#    #edge.attr.class[edge.attr.class=="character"]="char"
#    ## init node attributes
#    node.attr=list.vertex.attributes(refnet)
#    if (length(node.attr)){
#      node.attr.class = sapply(node.attr, class)
#      node.attr.class[node.attr.class=="character"]="char"
#      for (i in 1:length(node.attr))
#        net<- initNodeAttribute(net, attribute.name = node.attr[i],
#                              attribute.type = node.attr.class[i], default.value = "0")
#    }
#    ## our metagenomic network does not have edge attributes, set them all to 1
#    net <- initEdgeAttribute(net, attribute.name = "weight",
#                                attribute.type = "numeric", default.value = "1")
#    ## create a network window in Cytoscape
#    cw <- new.CytoscapeWindow ('net', graph = net, overwriteWindow = TRUE)
#    ## transmits the CytoscapeWindowClass's graph data, from R to Cytoscape, nodes,
#    ## edges, node and edge attributes
#    displayGraph (cw)
#  }

## ----echo=FALSE-----------------------------------------------
sessionInfo()

## ----closeConnetions------------------------------------------
allCon <- showConnections()
socketCon <- as.integer(rownames(allCon)[allCon[, "class"] == "sockconn"])
sapply(socketCon, function(ii) close.connection(getConnection(ii)) )

