#' Correction of batch effect in merged dataset
#' @description  This function performs various type of batch-correction on given merged
#' dataset and output batch-corrected data as a list.
#' @param experiment It is the merged dataset obtained at the previous step (`merge_experiments`).
#' @param model A R formula mentioning the biological factor to take into account during the correction.
#' @param k No. of confounders for RUV methods and denotes the "number of factors of unwanted variation to remove".
#' If k is not provided it will then estimated from input data.
#' @param batch meta-data column which has information about batch-label
#' @param filter A string. Should be one of following string- 'symbol', 'ensembl_gene_id', or 'entrezgene_id'
#' depending on gene label for the given dataset.
#' @return A list of the batch-corrected experiment.
#' @details  This function contains a cocktail of various batch-correction methods.
#'
#' For Microarray, `batch_correction` function has following methods-
#'
#'a) limma- Default limma setting,
#'b) GFS - Gene Fuzzy Score method.
#'c) Robust Quantile Normalization- quantile normalization method from PreprocessCore package,
#'d) ComBat- Implemented from SVA package,
#'e) Q_ComBat - It is Quantile+ parametric adjustment of ComBat,
#'f) MNN - is default mnnCorrect function from batchelor package.
#'g) naiveRandRUV_HK - default naiveRandRUV using Human HK genes from RUVnormalize package,
#'h) Q_naiveRandRUV_HK - is Qunatile+ naiveRandRUV_HK,
#'i) naiveRandRUV_empi.controls- default naiveRandRUV using empirically derived HK genes from RUVnormalize package,
#'j) Q_naiveRandRUV_empi.controls - is Qunatile+ naiveRandRUV_empi.controls.
#'
#' For bulk-RNAseq, `batch_correction` function has following methods-
#'
#' a) limma- default limma setting,
#' b) ComBat- implemented from SVA package,
#' c) ComBat_seq - implemented from ComBat-Seq package,
#' d) MNN- default mnnCorrect function from batchelor package,
#' e) RUVs- default RUVs function from RUVSeq package,
#' g) scBatch- batch correction method implemented in scBatch package.

#'
#' @import sva
#' @import RUVnormalize
#' @import limma
#' @import WGCNA
#' @import preprocessCore
#' @import Biobase
#' @import RUVSeq
#' @import batchelor
#' @import SummarizedExperiment
#' @import scBatch
#' @import edgeR

#' @include gfs_function.R
#' @include ComBat_seq.R
#' @include helper_ComBat_seq.R
#' @include Metadatacreation.R

#' @name batch_correction

data(control.genes, envir=environment())
#' @export
batch_correction <-function(experiment, model,k=NULL, batch="batch",filter)
{UseMethod("batch_correction")
}

#' @export
batch_correction.ExpressionSet<-function(experiment, model,k=NULL,batch="batch",filter)
  {
    edata <- experiment@assayData$exprs
    pheno <- experiment@phenoData
    model.data<-model.frame(model, experiment@phenoData@data[all.vars(model)])
    mod=model.matrix(model, data=model.data)


data.limma <-  ExpressionSet(assayData = removeBatchEffect(edata, batch =  experiment[[batch]],
                                                           design = mod), phenoData = experiment@phenoData)

data.GFS <-  ExpressionSet(assayData = gfs(edata), phenoData = experiment@phenoData)


data.quantile <-  normalize.quantiles.robust(edata,copy=TRUE,weights=NULL,remove.extreme="variance",
                                             use.median=FALSE,use.log2=TRUE) # using the geometric mean
colnames( data.quantile) <- colnames(edata)
rownames( data.quantile ) <- rownames(edata)
data.quantile <-  ExpressionSet(assayData = data.quantile, phenoData = experiment@phenoData)


data.ComBat1 <-  ExpressionSet(assayData = ComBat(dat=edata, batch= experiment[[batch]], par.prior=TRUE, prior.plots=FALSE),
                               phenoData = experiment@phenoData) #parametric adjustment
data.ComBat2 <-  ExpressionSet(assayData = ComBat(dat=edata, batch= experiment[[batch]], par.prior=FALSE, mean.only=TRUE),
                               phenoData = experiment@phenoData) #non-parametric adjustment, mean-only version
data.Q_ComBat <-  ExpressionSet(assayData = ComBat(dat=exprs(data.quantile), batch= experiment[[batch]], par.prior=TRUE, prior.plots=FALSE),
                               phenoData = experiment@phenoData) #parametric adjustment

data.MNN <-  ExpressionSet(assayData = as.matrix(assay(mnnCorrect(edata, batch= experiment[[batch]]))), phenoData = experiment@phenoData)

SV <- if(is.null(k))  {
  num_sv     <- sva::num.sv(dat = edata,mod = mod, method = "be")
  num_sv_l     <- sva::num.sv(dat = edata,mod = mod, method = "leek")
  SV <- ifelse(num_sv < num_sv_l,num_sv,num_sv_l)
} else {
  SV=k
}

empi.controls <- empirical.controls(edata,mod,mod0=NULL,n.sv=SV,type="norm")
genes <-  as.character(row.names(edata))
control.genes <- SelectBCM::control.genes

control.genes <-     if(filter=='symbol'){
  control.genes <- as.array(control.genes$hgnc_symbol)
} else if(filter=='entrezgene_id'){
  control.genes<- as.array(control.genes$entrezgene_id)
} else {
  control.genes <- as.array(control.genes$ensembl_gene_id)
}
ctrl.genes <- intersect(control.genes,genes)


data.naiveRandRUV_HK <-  ExpressionSet(assayData = t(naiveRandRUV(t(edata),ctrl.genes, nu.coeff=0.0, k=SV, tol=1e-6)),
                                       phenoData = experiment@phenoData)


data.Q_naiveRandRUV_HK <-  ExpressionSet(assayData = t(naiveRandRUV(t(exprs(data.quantile)),ctrl.genes, nu.coeff=0, k=SV, tol=1e-6)),
                                         phenoData = experiment@phenoData)


data.naiveRandRUV_empi.controls <-  ExpressionSet(assayData = t(naiveRandRUV(t(edata),empi.controls, nu.coeff=0, k=SV, tol=1e-6)),
                                                  phenoData = experiment@phenoData)


data.Q_naiveRandRUV_empi.controls <-  ExpressionSet(assayData = t(naiveRandRUV(t(exprs(data.quantile)),empi.controls, nu.coeff=0.0, k=SV, tol=1e-6)),
                                                    phenoData = experiment@phenoData)


batchcorrected <- list()
batchcorrected$data.limma <- data.limma
batchcorrected$data.GFS <- data.GFS
batchcorrected$data.quantile <- data.quantile
batchcorrected$data.ComBat1 <- data.ComBat1
batchcorrected$data.ComBat2 <- data.ComBat2
batchcorrected$data.Q_ComBat<- data.Q_ComBat

batchcorrected$data.mnncorrect <- data.MNN

batchcorrected$data.Q_naiveRandRUV_HK <- data.Q_naiveRandRUV_HK
batchcorrected$data.naiveRandRUV_HK <- data.naiveRandRUV_HK
batchcorrected$data.naiveRandRUV_empi.controls <- data.naiveRandRUV_empi.controls
batchcorrected$data.Q_naiveRandRUV_empi.controls <- data.Q_naiveRandRUV_empi.controls


batchcorrected
}


#' @export
batch_correction.SummarizedExperiment <-function(experiment, model, k=NULL, batch="batch")
  {
  log<-experiment@assays@data %>% names %>% switch(log_counts=TRUE, counts=FALSE)
  edata <- as.matrix(assay(experiment))
  edata <- edata[rowSums(edata)>0,]
  model.data<-model.frame(model, experiment@colData[all.vars(model)])
  mod=model.matrix(model, data=model.data)

  data.limma <-  SummarizedExperiment(assays=list(removeBatchEffect(edata, batch= as.numeric(as.factor(experiment$batch)),
                                                             design = mod)),
                                      colData = experiment@colData,
                                      metadata = experiment@metadata)

  data.ComBat <-  SummarizedExperiment(assays=list(ComBat(edata,batch= as.numeric(as.factor(experiment$batch)))),
                                      colData = experiment@colData,
                                      metadata = experiment@metadata)


  data.ComBat_seq_null <-  SummarizedExperiment(assays=list(ComBat_seq(edata,batch= as.numeric(as.factor(experiment$batch)))),
                                                colData = experiment@colData,
                                                metadata = experiment@metadata)
  data.ComBat_Seq_full <-  SummarizedExperiment(assays=list(ComBat_seq(edata,group= as.numeric(as.factor(experiment$Disease)),batch= as.numeric(as.factor(experiment$batch)))),
                                                colData = experiment@colData,
                                                metadata = experiment@metadata)

  data.MNN <-  SummarizedExperiment(assays=list(mnnCorrect(edata, batch=as.numeric(as.factor(experiment$batch)), cos.norm.out = F)@assays@data$corrected),
                                    colData = experiment@colData,
                                    metadata = experiment@metadata)

  SV <- if(is.null(k))  {
    num_sv     <- sva::num.sv(dat = edata,mod = mod, method = "be")
    num_sv_l     <- sva::num.sv(dat = edata,mod = mod, method = "leek")
    SV <- ifelse(num_sv < num_sv_l,num_sv,num_sv_l)
  } else {
    SV=k
  }


  data.RUVs = SummarizedExperiment(assays =list(RUVs(edata, cIdx=seq_len(nrow(edata)), k=SV,
                                                      scIdx=model.data %>% expand.grid %>% apply(1,paste) %>% makeGroups, isLog=TRUE)$normalizedCounts),
                                    colData = experiment@colData,
                                    metadata = experiment@metadata)



  ccc <- QuantNorm(edata,as.numeric(as.factor(experiment$batch)),logdat=F,method='row/column',cor_method='pearson',max=5)

  data.scBatch <-scBatchCpp(c=edata,d=ccc,m=5,w=diag(ncol(edata)),max=1200,tol=1e-10,step=0.0001,derif=scBatch::derif,verbose=F)
  colnames(data.scBatch) = colnames(edata)
  rownames(data.scBatch) = rownames(edata)

  data.scBatch = SummarizedExperiment(assays =list(data.scBatch),
                                   colData = experiment@colData,
                                   metadata = experiment@metadata)

 batchcorrected <- list()
  batchcorrected$data.limma <- data.limma
  batchcorrected$data.ComBat <- data.ComBat
  batchcorrected$data.ComBat_Seq_full <- data.ComBat_Seq_full
  batchcorrected$data.ComBat_seq_null <- data.ComBat_seq_null
  batchcorrected$data.MNN <- data.MNN
  batchcorrected$data.RUVs <- data.RUVs
  batchcorrected$data.scBatch  <- data.scBatch
  batchcorrected
}

