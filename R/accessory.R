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


