#' Assessment of batch-correction methods for bulk RNAseq experiments
#' @description  The function performs various evaluations on batch-corrected data and output
#' performance list of each individual method. Since there are number of ways batch-correction can
#' be evaluated and each method has its own limitation, we have used a cocktail of methods to
#' perform analysis.
#' @param  result A list of wrapped batch-corrected experiments obtained from the last step ('batch_correction').
#' @param  batch.factors A list of factors to perform PVCA analysis. Along with the batch factor, one biological factor
#' which can be used to assess over-fitting should be provided.
#' @param experiment A merged experiment without batch correction obtained from step ('merge_experiments').
#' @param N1 is the number of randomly picked samples for the BatchEntropy function.
#' @param N2 is the number of nearest neighbors of a sample (from all batches) to check (for BatchEntropy function).
#' @param filter A string. Should be one of following string- 'symbol', 'ensembl_gene_id', or 'entrezgene_id'
#' depending on gene label for the given dataset.
#' @examples  batch_evaluation(result, batch.factors=c("batch","sex"),experiments,'ensembl_gene_id')
#' @return A list of the evaluation methods on the batch-corrected experiment.

#' @import lme4
#' @import statmod
#' @import tibble

#' @importFrom  purrr map map_dfr
#' @importFrom  magrittr extract
#' @importFrom RANN nn2
#' @importFrom data.table as.data.table

#'
#' @include entropy.R
#' @include accessory.R
#' @include PVCA.R
#' @include pcRegression.R
#' @include batch_sil.R

#' @name batch_evaluation.RNAseq


#' @export
batch_evaluation.RNAseq <-function(result, batch.factors, experiment,N1,N2,filter)
{

  microarray.data <- function(filename) {
    i <- as.matrix((filename@assayData$exprs))}

  seq.data <- function(filename) {
    i <- as.matrix((filename@assays@data[[1]]))
  }

  transpose_df <- function(df) {
    t_df <- data.table::transpose(df)
    colnames(t_df) <- rownames(df)
    rownames(t_df) <- colnames(df)
    t_df <- t_df %>%
      tibble::as_tibble(.)
    return(t_df)
  }


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
                 as.numeric(experiment[["batch"]]),L=100, M=N1, k=N2)} ## Default values for Entropy calculation


  expression.input <- if(experiment@class =="SummarizedExperiment") {
    map(result,  seq.data)
  } else {
    map(result, microarray.data)
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



  silhouette.data <- map(expression.input, Batch.silhouette) %>% map_dfr(extract) %>% set_rownames("silhouette")
  pcRegression.data <- map(expression.input, Batch.pcRegression) %>% map_dfr(extract) %>% set_rownames("pcRegression")
  entropy.data <- map(expression.input, entropy)  %>% map_dfr(colMedians) %>% set_rownames("Entropy")
  pvca.data <- map(expression.input, PVCA_result.eset)  %>% map_dfr(extract,batch.factors)%>% transpose_df  %>% purrr::set_names(colnames(silhouette.data)) %>% tbl_df %>% set_rownames(batch.factors)




  GS.hvg <- as.data.table(t(raw.input)) %>%
    split(as.factor(experiment[["batch"]])) %>% lapply(data.frame)  %>% map((hvg_result_batches))
  GS.hvg <- GS.hvg %>% map('HVG') %>% bind_rows
  rownames(GS.hvg) <- row.names(raw.input)

  GS1 <- GS.hvg %>%
    rownames_to_column('gene') %>%
    filter_if(is.numeric, all_vars(. > 0)) %>%
    column_to_rownames('gene')
  GS1 <- as.vector(row.names(GS1))
  GS2 <- GS.hvg %>%
    rownames_to_column('gene') %>%
    filter_if(is.numeric, any_vars(. == 1)) %>%
    column_to_rownames('gene')
  GS2 <- as.vector(row.names(GS2))

  GS_corrected <- map(expression.input,hvg_result_corrected) %>% map('HVG') %>% bind_rows
  GS_corrected$name <- row.names(raw.input)

  GS_intersection <- GS_corrected %>% filter(name %in% GS1) %>% select(-name) %>% t
  GS_union <- GS_corrected %>% filter(name %in% GS2) %>% select(-name) %>% t

  #HVG.intersection <- as.data.frame(rowSums(GS_intersection)/length(GS1)) %>% t %>% tbl_df  %>% set_rownames("HVG.intersection")
  HVG.union <-  as.data.frame(rowSums(GS_union)/length(GS2)) %>% t %>% tbl_df%>% set_rownames("HVG.union")


  evaluation <- list()
  evaluation$pvca <- pvca.data
  evaluation$silhouette <- silhouette.data
  evaluation$pcRegression <- pcRegression.data
  evaluation$entropy <-  entropy.data
  evaluation$HVG.union<-  HVG.union
  evaluation
}

