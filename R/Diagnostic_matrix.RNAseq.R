#' Ranking of batch-correction methods for bulk RNAseq experiments
#' @description The function performs ranking of the evaluation data obtained at
#' the previous step 'batch_evaluation’ for the bulk RNAseq experiments. Here, methods are ranked on their individual performance,
#' and finally 'sumRank' is the final Rank of each method for the given dataset.
#' Rank 1 will be the best performer method.

#' @param  evaluation evaluation is the evaluation list obtained from the previous step `evaluation_matrix`.

#' @examples  diagnostic.RNAseq(assessment)
#' @return A list of two data-frame- (a) raw - simple data-frame output of evaluation matrix and, (b) ranked-
#' Ranked data-frame of evaluation matrix which has additional column `sumRank` containing final Rank of each method.
#' Ranks are in descending performance order, i.e. method having score 1 will be the best method.

#' @name diagnostic_matrix.RNAseq

#' @import GGally
#' @import Cairo
#' @importFrom dplyr dense_rank
#' @importFrom dplyr mutate_at
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggarrange
#' @importFrom colorspace divergingx_hcl
#' @importFrom hrbrthemes theme_ipsum

#' @export
diagnostic_matrix.RNAseq <-function(evaluation)

{

  minus <- function(x, na.rm = FALSE) (-x)
  mycolors <- divergingx_hcl(9, palette = "Roma")


  w_m <- do.call(rbind, lapply(evaluation, as.matrix)) %>% t %>% as.data.frame()
  w3_m <- w_m %>% mutate_at(vars(-c(batch,pcRegression,silhouette)),minus)
  w3_m <- mutate_each(w3_m, list(dense_rank))

  w3_m <- w3_m %>% mutate(sumRank = rowSums(.)) %>% mutate(across(c("sumRank"), dense_rank))
  row.names(w3_m) <- rownames(w_m)
  w3_m$method <- c("limma","ComBat","MNN","RUVs","scBatch","ComBat_seq_full",
                            "ComBat_seq_null","uncorrected","uncorrected1")

Rank <- list()
Rank$raw <- w_m
Rank$raw_rank <- w3_m
Rank

}

