#' Ranking and plotting of assessment of batch-correction methods
#' @description The function performs ranking and plotting of the evaluation data obtained at
#' the previous step 'batch_evaluation’. Here, methods are ranked on their individual performance,
#' and finally 'sumRank' is the final Rank of each method for the given dataset.
#' Rank 1 will be the best performer method.

#' @param  evaluation evaluation is the evaluation list obtained from the previous step `evaluation_matrix`.

#' @examples  diagnostic.RNAseq(assessment)
#' @return A list of two data-frame- (a) raw - simple data-frame output of evaluation matrix and, (b) ranked-
#' Ranked data-frame of evaluation matrix which has additional column `sumRank` containing final Rank of each method.
#' Ranks are in descending performance order, i.e. method having score 1 will be the best method. This function also output
#' `diagnostic plot`, where the x-axis is the evaluation methods and the y-axis is the Rank of each normalization method.
#' @name diagnostic.RNAseq

#' @import GGally
#' @importFrom dplyr dense_rank
#' @importFrom dplyr mutate_at

#' @export
diagnostic.RNAseq <-function(evaluation)

{

  minus <- function(x, na.rm = FALSE) (-x)
  mycolors <- divergingx_hcl(9, palette = "Roma")


  w_m <- do.call(rbind, lapply(evaluation, as.matrix)) %>% t %>% as.data.frame()
  w3_m <- w_m %>% mutate_at(vars(-c(batch,pcRegression,silhouette)),minus)
  w3_m <- mutate_each(w3_m, list(dense_rank))

  w3_m <- w3_m %>% mutate(sumRank = rowSums(.)) %>% mutate(across(c("sumRank"), dense_rank))
  row.names(w3_m) <- rownames(w_m)
  w3_m$method <- row.names(w3_m)


  p_ <- GGally::print_if_interactive

  svg("Diagnostic_plot.svg")

  parallelCoordinatePlot <-
    ggparcoord(data=w3_m, columns=1:6, groupColumn = "method",
               scale = "globalminmax",
               showPoints = TRUE,
               alphaLines = 1,splineFactor=10
    )  + scale_y_reverse()+scale_color_manual(values=mycolors)+
    theme_ipsum()+
    xlab("Evaluation method") +
    ylab( "Rank" ) +
    theme(legend.position="top") + theme(text = element_text(family = "Times New Roman"))+
    theme(axis.text.x = element_text(size=10, face="bold", color = "black"),
          axis.text.y = element_text(size=10, face="bold", color = "black"),
          legend.title=element_text(size=14, face="bold", color = "black"))



  p_(parallelCoordinatePlot)

  dev.off()


Rank <- list()
Rank$raw <- w_m
Rank$raw_rank <- w3_m
Rank



}

