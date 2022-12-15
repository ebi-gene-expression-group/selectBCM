#' Loading of experiments in a list
#' @description  Loading data in R- If the data files are already on your computer, you can use this step.
#' @param directory A character string that is the path to the directory where all the experiments to load, and not any other file, are stored in Rdata format.
#' @param item.SimpleList If the experiments are wrapped in a SimpleList object, this is the accessor (character string or integer) of the actual experiment object within this SimpleList. The experiment is unwrapped from this SimpleList.
#' @return A list of the loaded experiments.
#' @details  The experiment should be either bulk-RNAseq(Summarized experiment) or microarray(ExpressionSet).
#'WARNING: Different types of experiment should not be processed together. For example, one Summarized experiment with expressionSet experiment, otherwise it will cause an error in further steps.
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph intersection
#' @import stringr
#' @import Biobase
#' @import purrr
#' @import magrittr
#' @import ggplot2
#' @import gtools
#' @import readr
#' @import Biobase
#' @import tibble
#' @importFrom magrittr extract
#' @importFrom utils download.file
#' @import SummarizedExperiment
#' @export
load_experiments <- function(directory, names=dir(directory), item.SimpleList='rnaseq')
{
  directory %>% dir %>%
    map(~get(load(paste0(directory,'/',.x)))) %>%
    map(~if(class(.x)=='SimpleList') .x[[item.SimpleList]] else .x) %>%
    set_names(names)
}

#' download_experiments_from_ExpressionAtlas
#' @description  To use data from Expression Atlas that can be downloaded in `.Rdata` format.
#' @param experiments Character list of Expression atlas experiment id
#' @examples
#' experiments <- download_experiments_from_ExpressionAtlas('E-ERAD-169','E-GEOD-73175','E-MTAB-2801')
#' @return A list of the downloaded experiments.
#' @details  The experiment should be either bulk-RNAseq(SummarizedExperiment) or microarray(ExpressionSet).
#'WARNING: Different types of experiment should not be processed together. For example, one SummarizedExperiment with expressionSet experiment, otherwise it will cause an error in further steps.
#' @export

download_experiments_from_ExpressionAtlas<-function(..., destdir=getwd() %>% paste('experiments',sep='/')){
  if(!(destdir %>% dir.exists) ){
    destdir %>% dir.create
  }
  for(experiment in list(...)){
    paste0('http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/',experiment,'/',experiment,'-atlasExperimentSummary.Rdata') %>% download.file(destfile = paste0(destdir,"/",experiment,".Rdata"))
  }
  destdir %>% load_experiments
}

#' Removal of isolated experiments
#' @description  To correct the batch effect, one needs to take the biological characteristics of the samples into account. If no sample of an experiment shares biological characteristics with samples from other batches, it is not possible to correct the batch effect with these batches since one cannot distinguish the biological difference from the artifact. The function remove_isolated_experiments removes the isolated experiments and plots graphs of intersections between the experiments before and after removal.
#' @param experiments A list of wrapped experiments.
#' @param  biological.group  A character string indicating the biological covariate that makes the biological groups. This must be the name of a column of the experiment.
#' @example experiments %<>% remove_isolated_experiments('organism_part')
#' @return A list of wrapped experiments.
#' @details The experiments of the list that have no biological group in common with the other experiments are removed. This is an essential step to the integration of experiments using batch effect correction. Without it, the biological groups that are unshared between the experiments are confounding factors with the batch effect, and the latter cannot be corrected. The ultimate condition for this is the graph of biological intersections between the experiments to be a connected graph.
#'
#' WARNING: this function only checks and removes the experiments that are isolated. However, there can remain several clusters of experiments that are part of a connected subgraph, while the graph of all the experiments isn't connected. The graph of intersections after the removal of experiments is displayed so that the user can notice it. In this condition, the user has to choose manually some experiments that make a connected subgraph.
#' @export
remove_isolated_experiments<-function(experiments, biological.group){
  warning('The following batches were removed as they do not have a ',biological.group,' column:\n',
          names(experiments)[experiments %>% map(~is.null(.x[[biological.group]])) %>% unlist])
  experiments[experiments %>% map(~is.null(.x[[biological.group]])) %>% unlist]<-NULL
  batch<-experiments %>% imap(~.y %>% rep(dim(.x)[2])) %>% unlist(use.names=FALSE)
  group<-experiments %>% map(~.[[biological.group]]) %>% unlist(use.names=FALSE)
  groups <- group %>% split(batch)
  intersections<-NULL
  for(i in groups){
    for(j in groups){
      intersections%<>%c(length(intersect(i,j)))
    }
  }
  intersections%<>%matrix(length(groups))%<>%set_colnames(names(groups))%<>%set_rownames(names(groups))
  intersections %>% graph_from_adjacency_matrix %>% plot
  experiments[rownames(intersections)[rowSums(intersections!=0)<=1]]<-NULL
  message('The following batches were removed as they do not share any common point on ',biological.group,' column:\n',
          rownames(intersections)[rowSums(intersections!=0)<=1])
  batch<-experiments %>% imap(~.y %>% rep(dim(.x)[2])) %>% unlist(use.names=FALSE)
  group<-experiments %>% map(~.[[biological.group]]) %>% unlist(use.names=FALSE)
  groups <- group %>% split(batch)
  intersections<-NULL
  for(i in groups){
    for(j in groups){
      intersections%<>%c(length(intersect(i,j)))
    }
  }
  intersections%<>%matrix(length(groups))%<>%set_colnames(names(groups))%<>%set_rownames(names(groups))
  intersections %>% graph_from_adjacency_matrix %>% plot
  return(experiments)
}

