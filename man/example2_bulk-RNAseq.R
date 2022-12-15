---
  title: "Batchevaluation wrapper short introduction"
author: "Madhulika Mishra, Guillaume HEGER"
date: "01/06/2020"
---

  # Introduction
  # Over time, a vast amount of genomic information for the given phenotype has been accumulated in the public repository mainly in GEO and ArrayExpress. However, it is still challenging to integrate different high throughput experiments reflecting similar phenotypes because of various other non-biological confounding factors such as type of array, date of experiment, a laboratory where data was generated, etc and can be summarised as ‘batch effects’. To ComBat this issue, various batch-effect correction algorithms(BECAs) such as ComBat, SVA, and RUV have been developed to remove such batch effects from the integrated data and have shown promising results in mining biological signals. Evaluation of batch correction protocols are mainly looking at the Principle component plots or RLE plot or sometimes by measuring batch entropy.
  # This is a project to increase the reusability and reproducibility of these workflows, to facilitate batch correction, comparison, and benchmarking for ExpressionAtlas. Our strategy is to provide command-line access to individual library functions through simple wrapper scripts packaged in R wrapper. Batchevaluation R Wrapper contains implements of a variety of methods for batch correction of microarray and bulk RNA-seq data.

# Project overview
  # Batchevaluation analyses include step such as:

  #   * Fetching experiments from ExpressionAtlas based on a biological factor of interest, such as Disease, species,and tissue
  # * Removing isolated experiments
  # * Merging experiments in a single dataset
# * Evaluation of batch effect in the merged dataset
# * Correcting batch effects
# * Evaluation of batch correction methods


# These steps may be implemented in a variety of ways including stand-alone tools, scripts, or R package functions. To compare equivalent logical steps between workflows, and to ‘mix and match’ those components for optimal workflows is, therefore, a challenging exercise without additional infrastructure. The current R package provides flexibale framework to perform batch correction and evaluate effect of batch correction on the given dataset.


## Installation

# Installation should take less than 5 min.

### Via Github and devtools

# If you want to install the package directly from Github, I recommend to use the `devtools` package.

```R
library(devtools)
install_github('XXXXXXXX') link is not up yet
```

### Manually

# Please download the package as `zip` archive and install it via

```R
install.packages('Batchevaluation.zip', repos = NULL, type = 'source')
```

## Installing the dependencies (first use)
# To run the following functions, you will need some packages that can be installed using these commands in R :
  ```r
install.packages(c('magrittr','stringr','purrr','ggplot2','igraph','gtools'))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('sva')
BiocManager::install('RUVSeq')
BiocManager::install('batchelor')
```
## Loading data in R
# If the data files are already on your computer, you can use this step. If you want to download data from Expression Atlas, skip this part and go directly to "Downloading data from Expression Atlas".

# - Put the data files in `.Rdata` format in a directory containing only them. (The R objects contained in those files must be of class `SimpleList` with a `$rnaseq` slot which is a `SummarizedExperiment` object containing the data.)
# - Load all the experiments in a list using the function `load_experiments`

```r
experiments <- load_experiments('directory_path')
experiments<-load_experiments('./')
```
# Example Microarray data is provided in **example1** folder of the package.

## Downloading data from Expression Atlas
# If you want to use data from Expression Atlas that can be downloaded in `.Rdata` format, you can use the function `download_experiments_from_ExpressionAtlas` in this way :

  ```r
experiments <- download_experiments_from_ExpressionAtlas('E-MTAB-3718','E-MTAB-3725','E-GEOD-44366','E-GEOD-74747')
```

# This downloads the experiments in a new directory called "experiments" in your working directory and loads all the experiments in R within a list, using `load_experiments` function.

# After having loaded the experiments, you get a list of `SummarizedExperiment` or list of `ExpressionSet` object :
  ```
> experiments
$GSE1297_rawdata.Rdata
ExpressionSet (storageMode: lockedEnvironment)
assayData: 8872 features, 29 samples
element names: exprs
protocolData: none
phenoData
rowNames: 1 2 ... 3 (29 total)
varLabels: ID Assay ... Batch (9 total)
varMetadata: labelDescription
featureData: none
experimentData: use 'experimentData(object)'
Annotation:

  $GSE1711_rawdata.Rdata
ExpressionSet (storageMode: lockedEnvironment)
assayData: 12386 features, 24 samples
element names: exprs
protocolData: none
phenoData
sampleNames: SAMPLE212704SUB4576 SAMPLE212705SUB4576 ... SAMPLE212727SUB4576 (24 total)
varLabels: AtlasAssayGroup age ... biosource_type (9 total)
varMetadata: labelDescription
featureData: none
experimentData: use 'experimentData(object)'
Annotation:
  ```

## Removing the isolated experiments
To correct batch effect, one needs to take the biological characteristics of the samples into account (organism part in our example). If no sample of an experiment shares biological characteristics with samples from other batches, it is not possible to correct batch effect with these batches since one cannot distinguish the biological difference from the artifact. The function `remove_isolated_experiments` removes the isolated experiments and plots graphs of intersections between the experiments before and after removal.

```r
experiments %<>% remove_isolated_experiments('Tissue')
```

# WARNING : this function only removes the isolated experiments. Although it is still possible that two or more unconnected groups of experiments remain, within which the experiments are connected. In this case, batch effect correction is not possible neither and one has to choose a group of experiments manually.

# The two following plots are displayed by the function. The first one shows the graph of intersections of all the experiments before the removal of isolated ones. The second shows the same graph after their removal.

## Merging experiments in a single dataset
The function `merge_experiments` merges all the experiments in the list in a single `SummarizedExperiment` or `ExpressionSet` object and doesn't perform any correction. This function has two additional arguments `log` and `filter` (respectively set to `TRUE` and `FALSE` by default).

*  The `log` argument determines whether to perform log transformation on the data (recommended).
*  The `filter` argument determines whether to filter genes for which all the samples of a batch have zero-counts. Set it to `TRUE` if you have issues in running ComBat at the next step.

```r
experiments %<>% merge_experiments
#OR
experiments %<>% merge_experiments(log=TRUE,filter=FALSE)
#or any other settings
```
`experiments` is now an only `ExpressionSet` object containing the information about batches both in its `@metadata` and `@colData` slots :
```
ExpressionSet (storageMode: lockedEnvironment)
assayData: 8872 features, 56 samples
  element names: exprs
protocolData: none
phenoData
  sampleNames: GSM21215.cel GSM21218.cel ... GSM697308.CEL (56 total)
  varLabels: ID Assay ... batch (10 total)
  varMetadata: labelDescription
featureData: none
experimentData: use 'experimentData(object)'
Annotation:
```

## Correcting batch effect
The function `batch_correction` perform various type of batch-correction on given merged datset and output batch-corrected data as list.

Short Deatil of methods implemented in `batch_correction` function is given below-

* For Microarray, `batch_correction` function has following batch correction methods-

    * **limma**- default limma setting,  more detail about method can be learn from Limma  [documentation](https://www.bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf).
    * **GFS**- Gene Fuzzy Score, more detail about method can be Learn from [publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1327-8).
    * **Robust quantile normalization**- quantile normalisation method from PreprocessCore package. More information can be obtained from their [documentation](https://www.bioconductor.org/packages//2.7/bioc/html/preprocessCore.html).
    * **ComBat**- implemented from SVA package. More information can be obtained from their [documentation](https://bioconductor.riken.jp/packages/3.0/bioc/html/sva.html).
In current `batch_correction` method, there are two versions of ComBat- 1) **ComBat1**- for parametric adjustment and, 2) **ComBat2** - for non-parametric adjustment, mean-only version.
    * **Q_ComBatt** - is Quantile+ parametric adjustment of ComBat.
    * **MNN** - default mnnCorrect function from batchelor package. More information can be obtained from their [documentation](https://bioconductor.org/packages/release/bioc/html/batchelor.html).
    * **naiveRandRUV_HK** - default naiveRandRUV using Human HK genes from RUVnormalize package. More information can be obtained from their [documentation](https://bioconductor.org/packages/release/bioc/html/RUVnormalize.html). Human Housekeeping gene list is from [publication](https://academic.oup.com/gigascience/article/8/9/giz106/5570567).
    * **Q_naiveRandRUV_HK** - is Qunatile+ naiveRandRUV_HK
    * **naiveRandRUV_empi.controls** - default naiveRandRUV using data-empirically derived HK genes from RUVnormalize package.
    * **Q_naiveRandRUV_empi.controls** -is Qunatile+ naiveRandRUV_empi.controls.


* For bulk-RNAseq, `batch_correction` function has following batch correction methods-

    * **limma**- default limma setting,  more detail about method can be learn from Limma  [documentation](https://www.bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf).
    * **ComBat**- implemented from SVA package. More information can be obtained from their [documentation](https://bioconductor.riken.jp/packages/3.0/bioc/html/sva.html).
    * **ComBat_seq** - implemented from ComBat-Seq package.It is RNAseq versio  of ComBat. More information can be obtained from this [publication](https://www.biorxiv.org/content/10.1101/2020.01.13.904730v1).
    * **MNN** - default mnnCorrect function from batchelor package. More information can be obtained from their [documentation](https://bioconductor.org/packages/release/bioc/html/batchelor.html).
    * **RUVs** - default RUVs function from RUVSeq package.Here, k=1 is used. More information can be obtained from their [documentation](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html).
    * **scBatch** - batch correction method implemented in scBatch package. More information can be obtained from their [publication](https://academic.oup.com/bioinformatics/article/34/15/2634/4916062).


This function has three arguments `experiment`, `model` and `batch`.

*  The `experiment` argument is the input merged dataset.
*  The `model` argument is a R formula mentioning the biological factor to take into account during correction.
* The `batch` argument is a meta-data column which has information about batch-label.

Since, this merged dataset had some issue while merging, phenodata was added sperately. Phenodata is also avilable with the package in `data` folder and can be accessed like this:

```r
load("~/Batchevaluation/data/pheno.example1.Rdata")
experiments <- ExpressionSet(exprs(experiments), phenoData=pheno.experiment1)

```

After assiging phenodata, batch correction can be performed like this:
```r
result <- ExpressionSet(exprs(experiments), phenoData=pheno.experiment1)
```
Result is the list of batch-corrected data.

```
> summary(result)
                                  Length Class         Mode
data.limma                        1      ExpressionSet S4
data.GFS                          1      ExpressionSet S4
data.quantile                     1      ExpressionSet S4
data.ComBat1                      1      ExpressionSet S4
data.ComBat2                      1      ExpressionSet S4
data.Q_ComBat                     1      ExpressionSet S4
data.mnncorrect                   1      ExpressionSet S4
data.Q_naiveRandRUV_HK            1      ExpressionSet S4
data.naiveRandRUV_HK              1      ExpressionSet S4
data.naiveRandRUV_empi.controls   1      ExpressionSet S4
data.Q_naiveRandRUV_empi.controls 1      ExpressionSet S4
```
## Evaluation of batch-correction


