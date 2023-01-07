#' Rank plot of batch-correction methods provided by SelectBCM tool
#' @import dplyr
#' @import ggplot2
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom colorspace divergingx_hcl

#' @param  diagnostic_matirx  diagnostic_matirx obtained from the previous step `diagnostic_matrix.microarray`.

#' @name bcm_ranking


#' @export
bcm_ranking <- function(rank_data){

  rank_data$raw_rank$sumRank <- as.numeric(rank_data$raw_rank$sumRank)
  x1 <- as.data.frame(rank_data$raw_rank)

  p_ <- GGally::print_if_interactive

  svg("Ranking_plot.svg")

   P1 <- ggplot(data=x1, aes(x=reorder(method, -sumRank), y = sumRank)) +
    geom_bar(stat="identity") + theme_minimal() + coord_flip()+
    xlab("Batch correction method") +
    ylab("Rank") + theme(text = element_text(family = "Times New Roman"))+
    theme(axis.text.x = element_text(size=10, face="bold", color = "black"),
          axis.text.y = element_text(size=10, face="bold", color = "black"))

   p_(P1)

  dev.off()

}
