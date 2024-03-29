% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation_RNAseq.R
\name{batch_evaluation.RNAseq}
\alias{batch_evaluation.RNAseq}
\title{Assessment of batch-correction methods for bulk RNAseq experiments}
\usage{
batch_evaluation.RNAseq(result, batch.factors, experiment, N1, N2, filter)
}
\arguments{
\item{result}{A list of wrapped batch-corrected experiments obtained from the last step ('batch_correction').}

\item{batch.factors}{A list of factors to perform PVCA analysis. Along with the batch factor, one biological factor
which can be used to assess over-fitting should be provided.}

\item{experiment}{A merged experiment without batch correction obtained from step ('merge_experiments').}

\item{N1}{is the number of randomly picked samples for the BatchEntropy function.}

\item{N2}{is the number of nearest neighbors of a sample (from all batches) to check (for BatchEntropy function).}

\item{filter}{A string. Should be one of following string- 'symbol', 'ensembl_gene_id', or 'entrezgene_id'
depending on gene label for the given dataset.}
}
\value{
A list of the evaluation methods on the batch-corrected experiment.
}
\description{
The function performs various evaluations on batch-corrected data and output
performance list of each individual method. Since there are number of ways batch-correction can
be evaluated and each method has its own limitation, we have used a cocktail of methods to
perform analysis.
}
\examples{
 batch_evaluation(result, batch.factors=c("batch","sex"),experiments,'ensembl_gene_id')
}
