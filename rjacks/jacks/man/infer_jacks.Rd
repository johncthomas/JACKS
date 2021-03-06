% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{infer_jacks}
\alias{infer_jacks}
\title{Infer JACKS decomposition of gRNA effect into a sample-specific gene effect, and sample-independent gRNA efficacy.}
\usage{
infer_jacks(expt, target_genes = NULL, n_iter = 5, reference_library = NA)
}
\arguments{
\item{expt}{\code{\link{SummarizedExperiment}} of replicate-averaged gRNA count log2-fold changes in treatment compared to control.}

\item{target_genes}{A vector of gene names to run JACKS on if not all genes are required. Default NULL (all genes will be used).}

\item{n_iter}{An integer number of iterations of JACKS inference performed. Default 5.}

\item{reference_library}{Name of the gRNA library to be used in inference ("avana","gecko2","yusa_v1o", or the path to a local grna results file). If this is specified, gRNA efficacies are not estimated, which greatly speeds up the computation. Default NULL (recalculate gRNA efficacies).}
}
\value{
\code{\link{SummarizedExperiment}} object endowed with JACKS output.
}
\description{
\code{infer_jacks()} takes a \link{SummarizedExperiment}, which must pass the \link{validate_logfchange_table} check,
and infers the decomposition of JACKS model for all desired genes. The experiment object is endowed with
JACKS output: the estimated gene essentialities, their uncertainties, and directional p-values
(nominal and FDR corrected) are stored in experiment metadata, the gRNA efficacies with uncertainties in
rowData, and prediction of the decomposition and its uncertainty as assays.
}
