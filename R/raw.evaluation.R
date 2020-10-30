#' Detecting batch-effect in raw merged dataset
#' @description This is an accessory function that performs a subset of evaluation tests of
#'  'evaluation_matrix' function and provides estimates whether the merged dataset obtained
#'  after 'merge_experiments' requires batch correction or not. A higher value of pvca.batch, silhouette,
#'  pcRegression, and entropy is indicative of batch-effects in a raw merged dataset without
#'  having any correction.

#' @param  result A merged experiment without batch correction obtained from step ('merge_experiments').
#' @param  batch.factors A list of factors to perform PVCA analysis. Along with the batch factor, one biological factor
#' which can be used to assess over-fitting should be provided.
#' @param experiment A merged experiment without batch correction obtained from step ('merge_experiments').
#' @param N1 is the number of randomly picked cells for the BatchEntropy function.
#' @param N2 is the number of nearest nearest neighbors of the sample (from all batches) to check (for BatchEntropy function).
#' @param filter A string. Should be one of following string- 'symbol', 'ensembl_gene_id', or 'entrezgene_id'
#' depending on gene label for the given dataset.
#' @examples  detect_effect(experiments,experiment = experiments, batch.factors=c("batch","Disease"),10,10,'symbol')
#' @return A list of the evaluation methods on the batch-corrected experiment.

#' @import lme4
#' @import statmod

#' @importFrom  purrr map map_dfr
#' @importFrom  magrittr extract
#' @importFrom RANN nn2
#' @importFrom  tibble rownames_to_column
#' @importFrom  tibble column_to_rownames
#' @importFrom  dplyr bind_rows
#' @importFrom  dplyr filter_if
#' @importFrom  dplyr all_vars
#' @importFrom  dplyr any_vars
#' @importFrom  dplyr filter
#' @importFrom  dplyr select
#' @importFrom  dplyr tbl_df
#'
#' @include entropy.R
#' @include accessory.R
#' @include PVCA.R
#' @include pcRegression.R
#' @include batch_sil.R

#' @name detect_effect

#' @export
detect_effect <-function(result,experiment, batch.factors, N1,N2,filter)
{
  R1 <- list()
  R1$raw <- result

  transpose_df <- function(df) {
    t_df <- data.table::transpose(df)
    colnames(t_df) <- rownames(df)
    rownames(t_df) <- colnames(df)
    t_df <- t_df %>%
      tibble::as_tibble(.)
    return(t_df)
  }


  microarray.data <- function(filename) {
    i <- as.matrix((filename@assayData$exprs))}

  seq.data <- function(filename) {
    i <- as.matrix((filename@assays@data[[1]]))
  }
  pca.data <- function(filename) {
    pca.data <- prcomp(transposed.exprs(filename) , center=TRUE)}
  transposed.exprs <- function(filename) {
    i <- as.matrix(t((filename)))}

  PVCA_result.eset <- function(filename) {
    pvca_re <- PVCA(as.matrix(filename), meta.eset, 0.6,FALSE)
  }


  Batch.silhouette <- function(filename) {
    experiment.silhouette <- batch_sil(pca.data(filename), as.factor(experiment[["batch"]]))}

  Batch.pcRegression <- function(filename) {
    batch.pca <- pcRegression(pca.data(filename),batch,n_top=20)
    result <- batch.pca$pcRegscale}

  entropy <- function(filename) {
    BatchEntropy(pca.data(filename)$x[,1:2],
                 as.factor(experiment[["batch"]]),L=100, M=N1, k=N2)} ## Default values for Entropy calculation

  sex.silhouette <- function(filename) {
    experiment.silhouette <- batch_sil(pca.data(filename), as.factor(experiment[["sex"]]))}

  expression.input <- if(experiment@class =="SummarizedExperiment") {
    map(R1,  seq.data)
  } else {
    map(R1, microarray.data)
  }

  raw.input <- if(experiment@class =="SummarizedExperiment") {
    as.matrix((experiment@assays@data[[1]]))
  } else {
    as.matrix((experiment@assayData$exprs))
  }

  pheno.input <- if(experiment@class =="SummarizedExperiment") {
    as.matrix(experiment@colData)
  } else {
    as.matrix(experiment@phenoData@data)
  }

  meta.eset <- pheno.input[,batch.factors]
  batch <- as.factor(meta.eset[,"batch"])


  pvca.data <- map(expression.input, PVCA_result.eset) %>% map_dfr(extract,batch.factors) %>%  transpose_df %>% set_rownames(batch.factors) %>% purrr::set_names(c("raw"))
  silhouette.data <- map(expression.input, Batch.silhouette) %>% map_dfr(extract) %>% set_rownames("silhouette")
  pcRegression.data <- map(expression.input, Batch.pcRegression) %>% map_dfr(extract) %>% set_rownames("pcRegression")
  entropy.data <- map(expression.input, entropy)  %>% map_dfr(colMedians) %>% set_rownames("Entropy")

  raw <- list()
  raw$pvca <- pvca.data
  raw$silhouette <- silhouette.data
  raw$pcRegression <- pcRegression.data
  raw$entropy <-  entropy.data
  raw
}

