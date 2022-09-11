# SelectBCM wrapper: Short introduction
"Madhulika Mishra, Lucas Barck, Pablo Moreno, Guillaume Heger, Janet M. Thornton, Irene Papatheodorou"


# Introduction
Over time, a vast amount of genomic information has been accumulated in the public repository mainly in GEO and ArrayExpress for the given phenotype. However, it is still challenging to integrate different high throughput experiments reflecting similar phenotype because of various other non-biological confounding factors such as type of array, date of experiment, the laboratory where data was generated, etc. These variations can be summarized as `batch effects`. In order to solve this issue, various batch-effect correction algorithms(BECAs) such as ComBat, SVA, and RUV have been developed to remove such batch effects from the integrated data and have shown promising results in mining biological signals. Evaluation of batch correction protocols involves mainly looking at the Principal Component Analysis (PCA) plot or Relative log expression (RLE) plot or sometimes by measuring Batch entropy.

# Project overview
SelectBCM analyses include step such as:

* Fetching experiments from ExpressionAtlas based on a biological factor of interest, such as disease, species and tissue
* Removing isolated experiments
* Merging experiments in a single dataset
* detetction of batch effect in the merged dataset
* Correcting batch effects
* Evaluation of batch correction methods


These steps may be implemented in a variety of ways including stand-alone tools, scripts or R package functions. To compare equivalent logical steps between workflows, and to ‘mix and match’ those components for optimal workflows is, therefore a challenging exercise without having additional infrastructure. The current R package provides a flexibility to perform batch correction and evaluation of various batch-correction methods for given dataset. Finally, it provide **performance rank** of each BECAs for the given dataset.


## Via Github and devtools

If you want to install the package directly from Github, I recommend using the `devtools` package. Package uses r-base=3.6.3, therefore it is advised to initiate conda environment with r-base=3.6.3. Wrapper uses "foreign"  package as dependency which is no more supported by cran, therefore it should be installed via conda before hand. 

```
conda create -p ./venv  -c conda-forge r-base=3.6.3
conda install -c conda-forge r-foreign
```
```R
library(devtools)
install_github('https://github.com/ebi-gene-expression-group/selectBCM')
```

##  Manually

Please download the package as `zip` archive and install it via

```R
install.packages('SelectBCM.zip', repos = NULL, type = 'source')
```
# Overview of steps available in SelectBCM
SelectBCM package has several steps ranging from meta-experiment creation to batch-effect evaluation step (`Figure1`). In the current wrapper, scripts are written in a user-friendly way. Short description of each step and example is given below -

![Workflow of the SelectBCM package.](/data/figure1.png)

## loading library
Sometime loading of 'magrittr','purrr' and 'dplyr' with `SelectBCM` package  is deprecated, therefore, it is recommended to load all of these together.

**Recommendation:** Computation of scBatch requires high memory allocation for Rcpp, therefore it is advised to increase R memory.

```r
library(magrittr)
library(purrr)
library(dplyr)
library(SelectBCM)
```

## Loading data in R
If the data files are already on your computer, you can use this step. If you want to download data from Expression Atlas, skip this part and go directly to "Downloading data from Expression Atlas".

**Recommendation:** Microarray input data should be preprocessed- appropriately back-ground corrected without any normalization step, probe-to-gene level mapped and log-transformed. It is also recommended to remove low-expressed genes from data if possible. For bulk-RNaseq data, it should be count matrix .


**Steps:**

- Put the data files in `.Rdata` format in a directory containing only them. (The R objects contained in those files must be either `SummarizedExperiment` or `ExpressionSet` experiment.)
- Load all the experiments in a list using the function `load_experiments` 

```r
experiments <- load_experiments('directory_path')
experiments<-load_experiments('./')
```
Example Microarray data is provided in **example1** folder of the package.

## Downloading data from Expression Atlas
If you want to use data from Expression Atlas that can be downloaded in `.Rdata` format, you can use the function `download_experiments_from_ExpressionAtlas` in this way :

```r
experiments <- download_experiments_from_ExpressionAtlas(E-MTAB-8549','E-MTAB-5060')
*** will download RNAseq experiments from expressionatlas.
```

This downloads the experiments in a new directory called **"experiments"** in your working directory and loads all the experiments in R within a list, using `load_experiments` function.
After having loaded the experiments, you will get a list of either `SummarizedExperiment` or list of `ExpressionSet` objects.

**Caution**: Avoid mixing experiments of `SummarizedExperiment` with `ExpressionSet`. Experiments can only belong from any one of the classes only. All the selected experiments should have same gene id format.


## Removing the isolated experiments
To correct the `batch effect`, one needs to take the biological characteristics of the samples into account (`Tissue` in our example). If no sample of an experiment shares biological characteristics with samples from other batches, it is not possible to correct batch effect with these batches because one cannot distinguish the biological difference from the artifact. The function `remove_isolated_experiments` removes the isolated experiments and plots graphs of intersections between the experiments before and after removal.

```r
experiments %<>% remove_isolated_experiments('Tissue')
```

**WARNING:** this function only removes the isolated experiment. Although it is still possible that two or more unconnected groups of experiments remain, within which the experiments are connected. In this case, batch effect correction is not possible and one will have to choose a group of experiments manually.

The two following plots are displayed by the function. The first one shows the graph of intersections of all the experiments before the removal of isolated ones. The second shows the same graph after their removal.

## Merging experiments in a single dataset
The function `merge_experiments` merges all the experiments in the list in a single `SummarizedExperiment` or `ExpressionSet` object and doesn't perform any correction. This function has two additional arguments `log` and `filter` (respectively set to `TRUE` and `FALSE` by default).

*  The `log` argument determines whether to perform log transformation on the data (recommended for bulk RNAseq).
*  The `filter` argument determines whether to filter genes for which all the samples of a batch have zero-counts. Set it to `TRUE` if you have issues in running ComBat at the next step.

```r
experiments %<>% merge_experiments
#OR
experiments %<>% merge_experiments(log=TRUE,filter=FALSE)

```
`experiments` is now an only `ExpressionSet` object containing the information about batches both in its `@metadata` and `@colData` slots :

**Caution:** Sometimes during merging experiment, phenodata (SDRF) file gets corrupted, hence, it is advised to always check meta-data before proceeding further.


## Correcting batch effect
The function `batch_correction` performs various methods of batch-correction on a given merged dataset and output batch-corrected data as a list.

Short detail of methods implemented in `batch_correction` function are given below-

* For Microarray, `batch_correction` function has following batch correction methods-

    * **limma**- default limma setting,  more detail about the method can be learnt from Limma  [documentation](https://www.bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf).
    * **GFS**- Gene Fuzzy Score, more detail about the method can be learnt from [publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1327-8).
    * **Robust quantile normalization**- quantile normalisation method from PreprocessCore package. More information can be obtained from their [documentation](https://www.bioconductor.org/packages//2.7/bioc/html/preprocessCore.html).
    * **ComBat**- implemented from SVA package. More information can be obtained from their [documentation](https://bioconductor.riken.jp/packages/3.0/bioc/html/sva.html). 
In the current `batch_correction` method, there are two versions of ComBat- 1) **ComBat1**- for parametric adjustment and, 2) **ComBat2** - for non-parametric adjustment, mean-only version.
    * **Q_ComBat** - is Quantile+ parametric adjustment of ComBat.
    * **MNN** - default mnnCorrect function from batchelor package. More information can be obtained from their [documentation](https://bioconductor.org/packages/release/bioc/html/batchelor.html).
    * **naiveRandRUV_HK** - default naiveRandRUV using Human HK genes from RUVnormalize package. More information can be obtained from their [documentation](https://bioconductor.org/packages/release/bioc/html/RUVnormalize.html). The Human Housekeeping gene list is from this [publication](https://academic.oup.com/gigascience/article/8/9/giz106/5570567). No. of confounders were estimated using both "leek" and "be" method and then least no. of estimated confounders were used as input for all variants of **naiveRandRUV** method. 
    * **Q_naiveRandRUV_HK** - is Qunatile+ naiveRandRUV_HK
    * **naiveRandRUV_empi.controls** - default naiveRandRUV using data-empirically derived HK genes from RUVnormalize package.
    * **Q_naiveRandRUV_empi.controls** -is Qunatile+ naiveRandRUV_empi.controls.

* For bulk-RNAseq, `batch_correction` function has the following batch correction methods-

    * **limma**- default limma setting,  more detail about the method can be learnt from Limma  [documentation](https://www.bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf).
    * **ComBat**- implemented from SVA package. More information can be obtained from their [documentation](https://bioconductor.riken.jp/packages/3.0/bioc/html/sva.html). 
    * **ComBat_seq** - RNAseq version of ComBat, now implemented in `SVA` package. More information can be obtained from this [publication](https://www.biorxiv.org/content/10.1101/2020.01.13.904730v1).
    * **MNN** - default mnnCorrect function from batchelor package. More information can be obtained from their [documentation](https://bioconductor.org/packages/release/bioc/html/batchelor.html).
    * **RUVs** - default RUVs function from RUVSeq package.Here, k=1 is used. More information can be obtained from their [documentation](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html).
    * **scBatch** - batch correction method implemented in scBatch package. More information can be obtained from their [publication](https://academic.oup.com/bioinformatics/article/34/15/2634/4916062).
    
    
This function has following arguments-

*  The `experiment` argument is the input merged dataset.
*  The `model` argument is a R formula mentioning the biological factor to take into account during correction.
* The `batch` argument is a meta-data column which has information about batch-label. 
* The `k` argument is a no. of confounders for RUV methods and denotes the "number of factors of unwanted variation to remove". If k is not provided it will then be estimated from input data.
* The `filter` argument specifies the gene label for the given dataset. It Should be one of the following string- 'symbol', 'ensembl_gene_id' or 'entrezgene_id' depending on the gene label for the given dataset.


Since, this merged dataset had some issues while merging, phenodata was added separately. Phenodata is also available with the package in `data` folder and can be accessed like this: 

```r
load("~/SelectBCM/data/pheno.example1.Rdata")
experiments <- ExpressionSet(exprs(experiments), phenoData=pheno.experiment1)

``` 

After assigning phenodata, batch correction can be performed like this:
```r
result <- batch_correction(experiment= experiments,model=~Disease, batch = "batch",filter='symbol')
``` 
Result is the list of batch-corrected data:
```r
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

## Assessment of batch-correction methods

The function `batch_evaluation` performs various evaluations on batch-corrected data and output performance list of each method. Since there are no. of ways batch-correction can be evaluated and each method has some limitation, we have used a cocktail of methods to perform analysis. This function has both PCA-inspired as well as biology
inspired qualitative assessment protocol for batch-correction.

![Overview of implemented assessment methods.](/data/figure2.png)

Short detail of methods implemented in `batch_evaluation` function is given below-

  * **pvca**- A function for principal variance component analysis. The function is written based on the 'pvcaBatchAssess' function of the PVCA R package and slightly changed to make it more efficient and flexible for both microarray and bulk RNAseq gene-expression data. From https://github.com/dleelab/pvca.
  * **silhouette**- Determine `batch effect` using the silhouette coefficient (adopted from scone) with default setting (nPcs=3). Taken from `kBET` [package]( https://github.com/theislab/kBET).
  * **pcRegression**- Determines `batch effect` by a linear model fit of principal components and a batch (categorical) variable with a default setting (n_top=20). Taken from `kBET` [package]( https://github.com/theislab/kBET).
 * **entropy**- Determines the batch effect by computing the entropy of mixing. It is a parameter to quantify the extent of the intermingling of cells from different batches. Taken from `MNN` [package]( https://github.com/MarioniLab/MNN2017).For calculation, first two Pcs are used as input. Since,depending on no. of samples in merged dataset, it is important to choose N1 and N2 for Batchentropy calculation, we have put both arguments as variables here.
 * **HVG.union**-  This analysis calculates fraction of conserved highly variable genes (HVG) according to Brennecke et al., 2013 [publication](https://www.nature.com/articles/nmeth.2645). This calculation also accounts for perseverance of biological heterogeneity by measuring ratio of conservation of HVG genes after batch correction. This function calculates ratio of number of conserved HVG after batch correction/ Union of HVG genes among all the batches. Therefore, this function is less stringent compared to `HVG.intersection` function.

This function has following arguments-

*  The `result` argument is a list of wrapped batch-corrected experiments obtained from the last step ('batch_correction').
*  The `batch.factors` is a list of factors to perform `PVCA` analysis. Along with `batch` as factor, one `biological factor` which can be used to assess over-fitting should be provided. Providing more than one `biological factor` here, will create issues while plotting the results in next step.
* The `experiment` argument is the input merged dataset.
* The `N1` is the number of randomly picked samples for the BatchEntropy function.
* The `N2` is the number of nearest neighbours of the sample (from all batches) to check (for BatchEntropy function).
* The `filter` argument specifies the gene label for the given dataset. It Should be one of following string- 'symbol', 'ensembl_gene_id' or 'entrezgene_id' depending on the gene label for the given dataset.


Assessment of batch-corrected data can be performed like this:
```r
assessment <- batch_evaluation(result, batch.factors=c("batch","sex"), experiments,10,10,filter='symbol')
``` 
assessment is a nested list of evaluation scores for each of the evaluation protocols.


### Interpretation of result of evaluation protocols
The following table explains how results from each evaluation protocol should be interpreted. PVCA, silhouette index and PcRegression measures residual `batch-effect` in corrected data, therefore for these measurements lower the score, better the performance will be. HVG (HVG.union & HVG.inersection) measures inherent biological heterogeneity, therefore higher the score, better will be method. Entropy measures entropy of batch-mixing, therefore; higher the score better the method is.


## Diagnostic plot
Once, assessment is performed, in next step, results obtained from the `batch_evaluation` step can be further analyszed using
`Rank.plot` function which performs ranking and plotting of evaluation data obtained at the previous step. Here, methods are ranked on their performance and finally `sumRank` is the final Rank of each method for the given input dataset. Rank 1 will be the best performer method.

This function has only one argument- 

*  `diagnostic` is an evaluation list obtained from the previous step `batch_evaluation`.


```r
 final <- bcm_rank(assessment)
``` 
**final** is a list of two data-frame- (a) raw - simple data-frame output of evaluation matrix and, (b) ranked- Ranked data-frame of evaluation matrix which has additional column `sumRank` containing final Rank of each method. Ranks are in descending performance order, i.e. the method having score 1 will be the best method. This function also outputs a `diagnostic plot`, where x-axis is the evaluation protocol and y-axis is the Rank of each batch-correction method.

![Diagnostic plot.](/data/Diagnostic_plot.svg)



# Accessory function
## Detetction of batch-effects in raw merged dataset
This is an accessory function which performs a subset of evaluation tests of `batch_evaluation` function and provides estimates whether the merged dataset obtained after 'merge_experiments' requires batch correction or not. Higher value of pvca.batch, silhouette, pcRegression and entropy is a direct indicative of batch-effects in raw merged dataset.

This function has following arguments- 

*  The `result` argument is a is the merged dataset obtained from `merge_experiments` step.
*  The `batch.factors` is a list of factors to perform `PVCA` analysis. Along with `batch` as factor, one `biological factor` which can be used to assess over-fitting should be provided. Providing more than one `biological factor` here, will create issues while plotting the results.
* The `experiment` argument is again the same merged dataset.
* The `N1` is the number of randomly picked samples for the BatchEntropy function.
* The `N2` is the number of nearest neighbours of the sample (from all batches) to check (for BatchEntropy function).
* The `filter` argument specifies the gene label for the given dataset. It Should be one of the following string- 'symbol', 'ensembl_gene_id' or 'entrezgene_id' depending on the gene label for the given dataset.

Assessment of batch-effect on raw merged data can be performed like this:
```r
batch_effect.raw <- detect_effect(experiments,experiment = experiments, batch.factors=c("batch","Disease"),10,10,'symbol')
Result <- do.call(rbind, lapply(batch_effect.raw, data.frame))
``` 
Result:
```r
> Result
                   raw
pvca.batch   0.9994122
pvca.Disease 0.0000000
silhouette   0.7406456
pcRegression 0.9173469
entropy      0.0000000
``` 
