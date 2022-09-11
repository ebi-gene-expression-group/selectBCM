microarray.data <- function(filename) {
  i <- as.matrix((filename@assayData$exprs))}

seq.data <- function(filename) {
  i <- as.matrix((filename@assays@data[[1]]))
}

pca.data <- function(filename) {
  pca.data <- prcomp(transposed.exprs(filename) , center=TRUE)}
transposed.exprs <- function(filename) {
  i <- as.matrix(t((filename)))}


hvg_result_batches <- function(filename) {
  hvg <-   as.data.frame(find_hvg(dataframe=(as.matrix(t(filename))),  p.threshold=1e-2, plot=FALSE, return.ranks=FALSE,
                                  return.p=TRUE))
}

hvg_result_corrected <- function(filename) {
  hvg <-   as.data.frame(find_hvg(dataframe=(as.matrix((filename))),  p.threshold=1e-2, plot=FALSE, return.ranks=FALSE,
                                  return.p=TRUE))
}

access_data<-function(experiment) UseMethod("access_data")
access_data.SummarizedExperiment<-function(experiment){
  return(assay(experiment))
}
access_data.ExpressionSet<-function(experiment){
  return(experiment@assayData$exprs)
}

access_pheno<-function(experiment) UseMethod("access_pheno")
access_pheno.SummarizedExperiment<-function(experiment){
  return(experiment@colData)
}
access_pheno.ExpressionSet<-function(experiment){
  return(experiment@phenoData@data)
}

access_meta<-function(experiment) UseMethod("access_meta")
access_meta.SummarizedExperiment<-function(experiment){
  return(experiment@metadata)
}
access_meta.ExpressionSet<-function(experiment){
  return(experiment@protocolData@data)
}


expressed.data <- function(filename) {
  data <- assay(filename)  %>% as.data.frame()  %>% filter(rowSums(across(where(is.numeric)))!=0)

  se1 <- SummarizedExperiment(data, colData = filename@colData)

  return (se1)
}
remove_NA <- function(filename) {
  data <- assay(filename) %>% na.omit()
  se1 <- SummarizedExperiment(data, colData = filename@colData)
  return (se1)
}

expressed.data <- function(filename) {
  data <- assay(filename)  %>% as.data.frame()  %>% filter(rowSums(across(where(is.numeric)))!=0)

  se1 <- SummarizedExperiment(data, colData = filename@colData)

  return (se1)
}
remove_NA <- function(filename) {
  data <- assay(filename) %>% na.omit()
  se1 <- SummarizedExperiment(data, colData = filename@colData)
  return (se1)
}

