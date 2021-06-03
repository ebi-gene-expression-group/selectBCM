#' Ranking and plotting of assessment of batch-correction methods
#' @description The function performs ranking and plotting of the evaluation data obtained at
#' the previous step 'batch_evaluation’. Here, methods are ranked on their individual performance,
#' and finally 'sumRank' is the final Rank of each method for the given dataset.
#' Rank 1 will be the best performer method.

#' @param  evaluation evaluation is the evaluation list obtained from the previous step `evaluation_matrix`.

#' @examples  Rank.plot(assessment)
#' @return A list of two data-frame- (a) raw - simple data-frame output of evaluation matrix and, (b) ranked-
#' Ranked data-frame of evaluation matrix which has additional column `sumRank` containing final Rank of each method.
#' Ranks are in descending performance order, i.e. method having score 1 will be the best method. This function also output
#' `diagnostic plot`, where the x-axis is the evaluation methods and the y-axis is the Rank of each normalization method.

#' @import GGally
#' @importFrom dplyr dense_rank
#' @importFrom dplyr mutate_at

#' @export
Rank.plot <-function(evaluation)
{
  minus <- function(x, na.rm = FALSE) (-x)

evaluation <- if (evaluation$gender=="Gender-based silhoutte analysis is not performed as meta-experiment do not have a \"sex\" column"){
  evaluation[names(evaluation) != "gender"]
  } else
  {evaluation}

w_m <- do.call(rbind, lapply(evaluation, data.frame)) %>% t %>% as.data.frame

w3_m <- w_m %>% mutate_at(vars(-c(pvca.batch,pcRegression,silhouette)),minus)
w3_m <-across(w3_m, list(dense_rank))

w3_m <- w3_m %>% mutate(sumRank = rowSums(.)) %>% mutate_at(c("sumRank"), dense_rank)
row.names(w3_m) <- rownames(w_m)
w3_m$method <- row.names(w3_m)

p_ <- GGally::print_if_interactive

parallelCoordinatePlot <- if (c("gender") %in% names(data.frame(w3_m))=="TRUE"){
  ggparcoord( data=w3_m, columns=1:9,groupColumn = "method", scale = "globalminmax",
              showPoints=TRUE, alphaLines=0.4, splineFactor=10 ) + scale_y_reverse()+
    ggtitle( paste( "Diagnostic plot", sep="" ) ) +
    xlab("Evaluation method") +
    ylab( "Rank" ) +
    theme(legend.position="top")
} else
{ggparcoord( data=w3_m, columns=1:8,groupColumn = "method",scale = "globalminmax",
             showPoints=TRUE, alphaLines=0.4, splineFactor=10 ) + scale_y_reverse() +
    ggtitle( paste( "Diagnostic plot", sep="" ) ) +
    xlab("Evaluation method") +
    ylab( "Rank" ) +
    theme(legend.position="top") }

p_(parallelCoordinatePlot)

Rank <- list()
Rank$raw <- w_m
Rank$raw_rank <- w3_m
Rank
}

