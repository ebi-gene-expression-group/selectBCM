#' plotting of ranks obatined form diagnostic_matrix.microarray step
#' @description The function performs plotting of the evaluation data obtained at
#' the previous step 'diagnostic_matrix.microarray’.

#' @param  diagnostic_matirx  diagnostic_matirx obtained from the previous step `diagnostic_matrix.microarray`.

#' @examples  diagnostic_plot.microarray(diagnostic_matirx)
#' @return A list of two plots-
#'  (a) Output is `diagnostic_plot`, where the x-axis is the evaluation methods and the y-axis is the Rank of each normalization method.
#'   (b) Additional plot to summarise rankings for each method as a violin plot including summary statistics for each BCM.
#' @name diagnostic_plot.microarray

#' @import GGally
#' @import Cairo
#' @importFrom dplyr dense_rank
#' @importFrom dplyr mutate_at
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggarrange
#' @importFrom colorspace divergingx_hcl
#' @importFrom hrbrthemes theme_ipsum

#' @export

diagnostic_plot.microarray <-function(diagnostic_matirx)
{

  mycolors <- divergingx_hcl(11, palette = "Roma")

  w3_m <- diagnostic_matirx$raw_rank

  x11 <- as.data.frame((w3_m [1:6]))
  row.names(x11) <- w3_m$method
  x11 <- as.data.frame(t(x11[1:6]))
  plot_df <- pivot_longer(x11, cols = 1:11)


  diagnostic_plot <-
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


  diagnostic_summary <- ggplot(plot_df, aes(x=name, y=value, fill=name)) +  geom_violin(trim=TRUE, fill="gray")+
    ylim(0,10) +
    geom_boxplot(width=0.1)+ scale_color_manual(values=mycolors)+
    theme_ipsum()+
    xlab("Method") +
    ylab( "Rank" ) +
    theme(legend.position="top") + theme(text = element_text(family = "Times New Roman"))+
    theme(axis.text.x = element_text(size=10, face="bold", color = "black", angle=90,hjust=1),
          axis.text.y = element_text(size=10, face="bold", color = "black"),
          legend.title=element_text(size=14, face="bold", color = "black"))


  diagnostic_plot_list <<- list()

  diagnostic_plot_list$diagnostic_plot <- diagnostic_plot
  diagnostic_plot_list$diagnostic_stat_summary <- diagnostic_summary


  diagnostic_plot_list

}

