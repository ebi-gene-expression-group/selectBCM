BCM_ranking_microarray <- function(evaluation)
{UseMethod("BCM_ranking")
}

#' @export
  BCM_ranking_microarray.diagnostic <- function(evaluation)
{
  minus <- function(x, na.rm = FALSE) (-x)

  mycolors <- divergingx_hcl(11, palette = "Roma") 
  evaluation <- evaluation
  
  w_m <- do.call(rbind, lapply(evaluation, as.matrix)) %>% t %>% as.data.frame()
  w3_m <- w_m %>% mutate_at(vars(-c(batch,pcRegression,silhouette)),minus)
  w3_m <- mutate_each(w3_m, list(dense_rank))
  
  w3_m <- w3_m %>% mutate(sumRank = rowSums(.)) %>% mutate(across(c("sumRank"), dense_rank))
  row.names(w3_m) <- rownames(w_m)
  w3_m$method <- row.names(w3_m)
  
  p_ <- GGally::print_if_interactive
  
  svg("Diagnostic_plot.svg")
  
  parallelCoordinatePlot <- 
    ggparcoord(data=rank_data$raw_rank, columns=1:6, groupColumn = "method", 
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

#' @export
BCM_ranking_microarray.ranking <- function(rank_data){
  
rank_data$raw_rank$sumRank <- as.numeric(rank_data$raw_rank$sumRank)
x1 <- as.data.frame(rank_data$raw_rank)

 svg("Ranking_plot.svg")

 ggplot(data=x1, aes(x=reorder(method, -sumRank), y = sumRank)) +
  geom_bar(stat="identity") + theme_minimal() + coord_flip()+
  xlab("Batch correction method") +
  ylab("Rank") + theme(text = element_text(family = "Times New Roman"))+
  theme(axis.text.x = element_text(size=10, face="bold", color = "black"),
        axis.text.y = element_text(size=10, face="bold", color = "black"))

dev.off() 

}

