#' Rough merger of microarray experiments
#' @name merge_experiments.ExpressionSet
#' @description  The function merge_experiments merges all the experiments in the list in a single experiment object and doesn't perform any correction. This function has two additional arguments log and filter (set to TRUE by default).
#' It is recommended to give microarray data only with back-ground correction and after filtering low-expressed genes
#' @param  experiments A list of wrapped experiments.
#' @param  filter.unexpressed.genes A logical indicating if the genes that are unexpressed across all the samples of a batch shall be removed. This is for the proper functioning of the current version of ComBat in correct_batch_effect.
#' @param log A logical ndicating of whether the data shall be log-transformed(recommended).
#' @param force A logical. If the experiments of the list are of different classes, the function is stopped, unless force=TRUE.
#' @example experiments %<>% merge_experiments(log=TRUE,filter=FALSE)
#' @return Single ExpressionSet
#' @note Caution: Sometime during merging experiments, phenodata (SDRF) file gets corrupted, it is advised to always check meta-data before proceeding further.


#' @export
merge_experiment.ExpressionSet <- function(experiments, filter.unexpressed.genes=TRUE, log, force=FALSE)
{
  if(experiments %>% map(class) %>% unlist %>% unique %>% length %>% is_greater_than(1) & !force) stop("The experiments must have the same class. Their classes are :\n", experiments %>% map(class) %>% unlist)

  genes<-experiments %>% map(rownames)
  shared.genes<-genes %>% purrr::reduce(intersect)
  unshared.genes<-genes %>% map(setdiff %>% partial(y=shared.genes))

  data <-experiments %>% map(~.x %>% access_data)%>%  map(~.x,magrittr::extract(~.x,shared.genes,))
  data <- data %>% map (~.x %>% as_tibble(rownames = NA) %>% rownames_to_column()) %>% reduce(full_join,by ="rowname")
  geneid <- data$rowname
  data <- select(data,-"rowname")
  row.names(data) <- geneid

  warning(unshared.genes %>% unlist %>% unique %>% length,' genes have been removed as they are not shared across the batches.')
  batch<-experiments %>% imap(~.y %>% rep(ncol(.x))) %>% unlist(use.names=FALSE) %>% factor
  if(filter.unexpressed.genes){
    unexpressed.genes <- data %>% t %>% data.frame %>% split(batch) %>% map(~colSums(.)==0) %>% purrr::reduce(`&`)
    data <- data %>% filter(rowSums(across(where(is.numeric)))!=0)
    message(sum(unexpressed.genes),' genes have been removed as they were unexpressed across the samples of a batch.')
  }

  if(log) data %<>% log1p

  expressionset <- ExpressionSet(
    assayData=as.matrix(data),
    phenoData=experiments %>% map(~.x@phenoData@data) %>% smartbind(list=.) %>% set_rownames(experiments %>% map(colnames) %>% unlist) %>% cbind(batch) %>% AnnotatedDataFrame)

  return(expressionset)
}
