% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_functions.R
\name{cca_cv_boot}
\alias{cca_cv_boot}
\title{Cross-Validtation Bootstrap CCA}
\usage{
cca_cv_boot(
  X_FIT,
  Y_FIT,
  ncomp = 10,
  Nboot = 30,
  Nfolds = 10,
  ProcrustX = NULL,
  ProcrustY = NULL,
  UseProgressBar = TRUE,
  UseProcrustes = TRUE
)
}
\arguments{
\item{X_FIT}{Numeric Matrix or Data Frame [N, P1] containing the predictor variables.}

\item{Y_FIT}{Numeric Matrix or Data Frame [N, P2] containing the outcome variables.}

\item{ncomp}{Numeric Scalar. Number of CCA components to keep in analyses. Must be equal to or less than min(P1,P2).}

\item{Nboot}{Numeric Scalar. Number of times to repeat k-fold cross-validation (yes its confusing it says "boot").}

\item{Nfolds}{Numeric Scalar. Number of}

\item{ProcrustX}{Numeric Matrix [ncomp, P1] containing target matrix for Procrustes Analysis. All CCA predictor raw coefficients obtained during the bootstrap resampling will be rotated to this target matrix.}

\item{ProcrustY}{Numeric Matrix [ncomp, P2] containing target matrix for Procrustes Analysis. All CCA outcome raw coefficients obtained during the bootstrap resampling will be rotated to this target matrix.}

\item{UseProgressBar}{Logical. Whether to show progress bar.}

\item{UseProcrustes}{Logical. Whether to use procrustes analysis.}
}
\description{
Algorithm for using cross-validation & bootstrap resampling to find unbiased canonical correlation and their
sampling error. On each iteration, N-fold (default is 10 fold) cross-validation is used to generated predicted canonical
variates for the complete sample. Following this, the predicted variates are bootstrap resampled and
canonical correlations are estimated from them.
Because bootstrap resampling can change the order of canonical variates that are extracted, or sign
flipping can occur in some cases (i.e. a very similar latent variable is extracted but on some occasions
the loadings are mostly positive or negative), we rotate the loadings in each during cross-validation to
map onto the loadings generated from a smaller dataset (ProcrustX & ProcrustY)
}
