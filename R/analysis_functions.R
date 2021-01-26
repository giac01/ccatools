#' Canonical Correlations
#'
#' @description
#' Sightly quicker than stats::cancor, and allows you to easily fit cca models in one dataset and find predicted variates/canonical correlations in another dataset.
#'
#' The number of predictor (X) and outcome (Y) variables are denoted by P1 and P2 below, and N is sample size.
#'
#' The function also allows the user to rotate the raw coefficients using Procrustes Analyses to target matrices (ProcrustX & ProcrustY), prior to estimating variates.
#'
#' @param X_FIT Numeric Matrix [N, P1] containing the training dataset predictor variables.
#' @param Y_FIT Numeric Matrix [N, P2] containing the training dataset outcome variables.
#' @param X_PRED Numeric Matrix [N, P1] containing the testing dataset predictor variables. Variables should be ordered in the same way as for X_FIT.
#' @param Y_PRED Numeric Matrix [N, P2] containing the testing dataset outcome variables. Variables should be ordered in the same way as for Y_FIT.
#' @param ncomp Numeric Scalar. Number of CCA components to keep in analyses. Must be equal to or less than min(P1,P2).
#' @param ProcrustX Numeric Matrix [ncomp, P1] containing target matrix for Procrustes Analysis. All CCA predictor raw coefficients obtained during the bootstrap resampling will be rotated to this target matrix.
#' @param ProcrustY Numeric Matrix [ncomp, P2] containing target matrix for Procrustes Analysis. All CCA outcome raw coefficients obtained during the bootstrap resampling will be rotated to this target matrix.
#' @param SafetyChecks Checks the input provided for mistakes (default = FALSE).
#'
#' @return A list containing the following components
#' \itemize{
#'   \item xcoef - Estimated raw coefficients (CCA weights) for the x (predictor) variables.
#'   \item ycoef - Estimated raw coefficients (CCA weights) for the y (outcome) variables.
#'   \item variates - Variates (latent variable scores), estimated from the raw coefficient in X_FIT/Y_FIT.
#'   cbind.data.frame(xvariates, yvariates). If X_PRED and Y_PRED are provided, then  variates (and xvariates/yvariates) will be the predicted latent variable scores  from X_PRED/Y_PRED matrices.
#'   \item xvariates - Variates (latent variable scores) for predictor variables.
#'   \item yvariates - Variates (latent variable scores) for outcome variables.
#'   \item cc_pred - Predicted Canonical Correlations, estimated using from X_PRED & Y_PRED if provided. If no PRED matrices are specified, then this simply returns the estimated canonical correlations from the training datsets (X_FIT/Y_FIT).
#'   \item cc_fit - Estimated Canonical Correlation estimated from the training datasets (X_FIT/Y_FIT).

#' }
#'
#' @export
#'
.cca = function(X_FIT,Y_FIT,X_PRED=NULL,Y_PRED=NULL,
                  ncomp=10,
                  ProcrustX = NULL, ProcrustY = NULL,
                  SafetyChecks=TRUE){
  # browser()
  #Check some basic things
  if (SafetyChecks){
    if (!is.matrix(X_FIT) | !is.matrix(Y_FIT)) stop("both X_FIT and Y_FIT should be matrices")
    if (nrow(X_FIT)!=nrow(Y_FIT)) stop("nrow of X_FIT and Y_FIT do not match")
    if (!is.null(ProcrustX)){
      if (ncol(ProcrustX)!=ncomp) stop("ProcrustX should have same number of columns to ncomp")
    }
    if (!is.null(ProcrustY)){
      if (ncol(ProcrustY)!=ncomp) stop("ProcrustY should have same number of columns to ncomp")
    }
    if (!is.null(ncomp)){
      if ((ncomp>ncol(X_FIT)) | (ncomp>ncol(Y_FIT)) ) stop("ncomp should be equal to or less than smallest number of variables in X_FIT or X_FIT")
    }
  }

  #Set ncomp to the minimum of ncol
  ncomp = min(c(ncomp,ncol(X_FIT),ncol(Y_FIT)))

  #Scale FIT datasets
  X = as.matrix(Rfast::standardise(X_FIT, center = TRUE, scale = TRUE))
  Y = as.matrix(Rfast::standardise(Y_FIT, center = TRUE, scale = TRUE))

  #Scale the PRED matrix using the FIT m and sd  - X variables
  if (is.null(X_PRED)){
    X_PRED = X
  } else {
    # Scale the PRED matrix by the FIT matrix
    X_FIT_mat = as.matrix(X_FIT)
    X_Mean = matrixStats::colMeans2(X_FIT_mat)
    X_SD = matrixStats::colSds(X_FIT_mat)
    X_PRED = t((t(X_PRED)-X_Mean)/X_SD)
  }

  #Scale the PRED matrix using the FIT m and sd  - Y variables
  if (is.null(Y_PRED)){
    Y_PRED = Y
  } else {
    # Scale the PRED matrix by the FIT matrix
    Y_FIT_mat = as.matrix(Y_FIT)
    Y_Mean = matrixStats::colMeans2(Y_FIT_mat)
    Y_SD = matrixStats::colSds(Y_FIT_mat)
    Y_PRED = t((t(Y_PRED)-Y_Mean)/Y_SD)
  }

  # Estimate Eigenvalues and Eigenvectors...
  Sxx  = crossprod(X,X)
  Syx  = crossprod(Y,X)
  Syy  = crossprod(Y,Y)
  Sxy  = crossprod(X,Y)

  #Slightly slower to use below code!
  # Sxx  = t(X) %*% X
  # Syx  = t(Y) %*% X
  # Syy  = t(Y) %*% Y
  # Sxy  = t(X) %*% Y

  V = solve(Syy) %*% Syx %*% solve(Sxx) %*% Sxy
  W = solve(Sxx) %*% Sxy %*% solve(Syy) %*% Syx

  Eig_V = eigen(V)
  Eig_W = eigen(W)

  cc = sqrt(Eig_V$values)[1:ncomp] # Canonical correlations from FITTED MATRIX

  # Get loadings
  ycoef = apply(as.matrix(Eig_V$vectors[,1:ncomp]),2,as.numeric)
  xcoef = apply(as.matrix(Eig_W$vectors[,1:ncomp]),2,as.numeric)

  # Optional Rotation of Loadings

  if (!is.null(ProcrustX)){
    xcoef = as.matrix(MCMCpack::procrustes(xcoef  , ProcrustX)$X.new)

  }
  if (!is.null(ProcrustY)){
    ycoef = as.matrix(MCMCpack::procrustes(ycoef  , ProcrustY)$X.new)
  }

  # Estimate Canonical Variates

  yvariates =  Y_PRED %*% ycoef
  colnames(yvariates) = paste0("Y", 1:ncomp)

  xvariates =  X_PRED %*% xcoef
  colnames(xvariates) = paste0("X", 1:ncomp)

  variates = cbind.data.frame(xvariates, yvariates)

  cc_pred = Rfast::corpairs(xvariates, yvariates)

  # Return Output

  out = list(
    #Loadings
    ycoef = ycoef,
    xcoef = xcoef,

    #variates
    variates = variates,
    xvariates = xvariates,
    yvariates = yvariates,

    #Canonical Correlations
    cc_fit = cc,                            #Canonical correlations estimated from the input data (fitted canonical correlations)
    cc_pred = cc_pred                       #Canonical correlations estimated from the prediction data (predicted canonical correlations)
  )


  return(out)

}

#' Split-Half CCA code
#'
#' Run CCA model in training dataset, and and validate performance in testing dataset.
#' The function estimates confidence intervals and p-values using standard inferential methods of Pearson's correlation coefficient.
#'
#' @param X_FIT Numeric Matrix or Data Frame [N, P1] containing the training dataset predictor variables.
#' @param Y_FIT Numeric Matrix or Data Frame [N, P2] containing the training dataset outcome variables.
#' @param X_PRED Numeric Matrix or Data Frame [N, P1] containing the testing dataset predictor variables. Variables should be ordered in the same way as for X_FIT.
#' @param Y_PRED Numeric Matrix or Data Frame [N, P1] containing the testing dataset outcome variables. Variables should be ordered in the same way as for Y_FIT.
#' @param ncomp Numeric Scalar. Number of CCA components to keep in analyses. Must be equal to or less than min(P1,P2).
#' @param alpha Numeric Scalar. Alpha level for estimating a 100(1-alpha)\% confidence interval for each canonical correlation. Default is .05 for estimating a 95\% confidence interval.
#'
#' @export
#'
cca_splithalf = function(X_FIT,Y_FIT,X_PRED,Y_PRED,
                            ncomp=NULL, alpha = 0.05){

  if (is.data.frame(X_FIT)) X_FIT  = as.matrix(X_FIT)
  if (is.data.frame(Y_FIT)) Y_FIT  = as.matrix(Y_FIT)
  if (is.data.frame(X_PRED)) X_PRED = as.matrix(X_PRED)
  if (is.data.frame(Y_PRED)) Y_PRED = as.matrix(Y_PRED)

  model_results = .cca(X_FIT=X_FIT,Y_FIT=Y_FIT,X_PRED=X_PRED,Y_PRED=Y_PRED,
                         ncomp=ncomp,
                         ProcrustX = NULL, ProcrustY = NULL,
                         SafetyChecks=TRUE)

  # Estimate Confidence Intervals

  samplesize = nrow(X_PRED)
  predicted_cc = model_results$cc_pred
  confint_cc = sapply(predicted_cc, function(x) CorrelationCIEstimator(x, samplesize, alpha=alpha)$CI)
  rownames(confint_cc) = c("LB","UB") # 'lower bound' and 'upper bound' of confidence interval
  colnames(confint_cc) = paste0("cc",1:ncol(confint_cc))
  pvalue_cc = sapply(predicted_cc, function(x) CorrelationCIEstimator(x, samplesize, alpha=alpha)$p)
  combined_cc = rbind.data.frame(predicted_cc, confint_cc, pvalue_cc)
  rownames(combined_cc) = c("cc", "LB", "UB", "p")

  return(list(
    model_results = model_results,
    predicted_cc  = predicted_cc,
    confint_cc    = confint_cc,
    pvalue_cc     = pvalue_cc,
    combined_cc   = combined_cc
  ))

}

#' Percentile Bootstrap Estimation of canonical correlation coefficients
#'
#' This function runs the .cca Canonical Correlation Analysis function multiple times to assess variability in the CCA loadings and canonical correlations
#' Because bootstrap resampling can change the order of canonical variates that are extracted, or sign flipping can occur
#' in some cases (i.e. a very similar latent variable is extracted but on some occasions the loadings are mostly positive or negative), we rotate the loadings
#' in each bootstrap resample to map onto the loadings generated from the full, raw input datsets.
#'
#'
#' @param X_FIT Numeric Matrix or Data Frame [N, P1] containing the predictor variables.
#' @param Y_FIT Numeric Matrix or Data Frame [N, P2] containing the outcome variables.
#' @param ncomp Numeric Scalar. Number of CCA components to keep in analyses. Must be equal to or less than min(P1,P2).
#' @param Nboot Numeric Scaler. Number of times to repeat bootstrap resampling.
#' @param ProcrustX Numeric Matrix [ncomp, P1] containing target matrix for Procrustes Analysis. All CCA predictor raw coefficients obtained during the bootstrap resampling will be rotated to this target matrix.
#' @param ProcrustY Numeric Matrix [ncomp, P2] containing target matrix for Procrustes Analysis. All CCA outcome raw coefficients obtained during the bootstrap resampling will be rotated to this target matrix.
#'
#' @export
#'
coef_boot = function(X_FIT,Y_FIT, ncomp=10, Nboot=30,ProcrustX = NULL, ProcrustY = NULL){
  # browser()
  pb <- utils::txtProgressBar(min = 1, max = Nboot, style = 3)

  if (is.data.frame(X_FIT)) X_FIT  = as.matrix(X_FIT)
  if (is.data.frame(Y_FIT)) Y_FIT  = as.matrix(Y_FIT)

  #Run model on full dataset
  CCA_OriginalData = .cca(X_FIT=X_FIT, Y_FIT=Y_FIT, X_PRED=NULL, Y_PRED=NULL,
                            ProcrustX = ProcrustX, ProcrustY = ProcrustY,
                            ncomp=ncomp,
                            SafetyChecks = TRUE)

  #Lists to store bootstrap data in
  ycoef_ROTATED = list()  #Store rotated loadings from bootstrap resamples in here
  xcoef_ROTATED = list()
  cc = list()                #Store canonical correlations from bootstrap resamples in hre

  for(i in 1:Nboot){
    # Randomly bootstrap resample data and fit CCA model.
    BootResample = base::sample(nrow(X_FIT), replace = TRUE)


    #Run CCA on bootstrap resampled data, rotating loadings using procrustes to map  onto the original model
    CCA_BootData = .cca(X_FIT=X_FIT[BootResample,], Y_FIT=Y_FIT[BootResample,],
                          ProcrustX = CCA_OriginalData$xcoef, ProcrustY = CCA_OriginalData$ycoef,
                          X_PRED=NULL, Y_PRED=NULL,
                          ncomp=ncomp)

    xcoef_ROTATED[[i]] = CCA_BootData$xcoef
    ycoef_ROTATED[[i]] = CCA_BootData$ycoef
    cc[[i]]  = CCA_BootData$cc_pred

    utils::setTxtProgressBar(pb, value=i)

  }
  #Estimate quantiles from boostrap distribution

  xcoef_Quantiles = apply(base::simplify2array(xcoef_ROTATED), 1:2, stats::quantile, prob = c(0.025, .5, 0.975))
  xcoef_Quantiles = lapply(1:ncomp, function(i) xcoef_Quantiles[,,i])
  for(i in 1:ncomp){
    colnames(xcoef_Quantiles[[i]]) = colnames(X_FIT)
    xcoef_Quantiles[[i]] = data.frame(t(xcoef_Quantiles[[i]]))
    xcoef_Quantiles[[i]]$original = CCA_OriginalData$xcoef[,i]
  }

  ycoef_Quantiles = apply(base::simplify2array(ycoef_ROTATED), 1:2, stats::quantile, prob = c(0.025, .5, 0.975))
  ycoef_Quantiles = lapply(1:ncomp, function(i) ycoef_Quantiles[,,i])
  for(i in 1:ncomp){
    colnames(ycoef_Quantiles[[i]]) = colnames(Y_FIT)
    ycoef_Quantiles[[i]] = data.frame(t(ycoef_Quantiles[[i]]))
    ycoef_Quantiles[[i]]$original = CCA_OriginalData$ycoef[,i]
  }

  cc_quantiles = apply(base::simplify2array(cc), 1, stats::quantile, prob = c(0.025, .5, 0.975))


  # Ouput Data

  out = list(
    xcoef_Quantiles=xcoef_Quantiles,
    ycoef_Quantiles=ycoef_Quantiles,
    cc_quantiles,
    cc=cc
  )

  return(out)

}


#'  Cross-Validtation Bootstrap CCA
#'
#' Algorithm for using cross-validation & bootstrap resampling to find unbiased canonical correlation and their
#' sampling error. On each iteration, N-fold (default is 10 fold) cross-validation is used to generated predicted canonical
#' variates for the complete sample. Following this, the predicted variates are bootstrap resampled and
#' canonical correlations are estimated from them.
#' Because bootstrap resampling can change the order of canonical variates that are extracted, or sign
#' flipping can occur in some cases (i.e. a very similar latent variable is extracted but on some occasions
#' the loadings are mostly positive or negative), we rotate the loadings in each during cross-validation to
#' map onto the loadings generated from a smaller dataset (ProcrustX & ProcrustY)
#'
#' @param ProcrustX Numeric Matrix [ncomp, P1] containing target matrix for Procrustes Analysis. All CCA predictor raw coefficients obtained during the bootstrap resampling will be rotated to this target matrix.
#' @param ProcrustY Numeric Matrix [ncomp, P2] containing target matrix for Procrustes Analysis. All CCA outcome raw coefficients obtained during the bootstrap resampling will be rotated to this target matrix.
#' @param X_FIT Numeric Matrix or Data Frame [N, P1] containing the predictor variables.
#' @param Y_FIT Numeric Matrix or Data Frame [N, P2] containing the outcome variables.
#' @param ncomp Numeric Scalar. Number of CCA components to keep in analyses. Must be equal to or less than min(P1,P2).
#' @param Nboot Numeric Scalar. Number of times to repeat k-fold cross-validation (yes its confusing it says "boot").
#' @param Nfolds Numeric Scalar. Number of
#' @param UseProgressBar Logical. Whether to show progress bar.
#' @param UseProcrustes  Logical. Whether to use procrustes analysis.
#'
#' @export
#'
cca_cv_boot = function(X_FIT,Y_FIT, ncomp=10, Nboot=30, Nfolds=10,ProcrustX = NULL, ProcrustY = NULL, UseProgressBar=TRUE, UseProcrustes=TRUE){
  if (is.data.frame(X_FIT)) X_FIT  = as.matrix(X_FIT)
  if (is.data.frame(Y_FIT)) Y_FIT  = as.matrix(Y_FIT)

  # browser()
  if (UseProgressBar){
    pb <- utils::txtProgressBar(min = 0, max = Nboot, style = 3)
  }

  CCA_OriginalData = .cca(X_FIT=X_FIT, Y_FIT=Y_FIT, X_PRED=NULL, Y_PRED=NULL, ncomp=ncomp, ProcrustX = ProcrustX, ProcrustY = ProcrustY)

  if (UseProcrustes==FALSE){
    CCA_OriginalData$xcoef = NULL # By setting this to NULL, the .cca function called below will not use procrustes rotations
    CCA_OriginalData$ycoef = NULL
  }

  cc_CVboot = list()
  cc_CV = list()
  R2_matrix = list()

  for(b in 1:Nboot){
    #Divide data into folds...
    df_folds = caret::createFolds(1:nrow(X_FIT), k=Nfolds)
    # df_folds = cut(sample(nrow(X_FIT)), breaks=Nfolds, labels=FALSE)

    Variate_predictions = list()
    for (f in 1:Nfolds){

      # Fit_Index = (df_folds!=f) #Rows to select to fit data to
      # Pred_Index = (df_folds==f) #Rows to select to make predictions for
      #
      Fit_Index = base::unlist(df_folds[-f], FALSE, FALSE) #Row numbers - Training Data
      Pred_Index = base::unlist(df_folds[f], FALSE, FALSE) #Row numbers - Hold-Out Data

      #Estimate CCA in trainning dataset and generate list of predictions for hold-out data - and append to Variate_predictions list
      Variate_predictions =  c(Variate_predictions,
                               list(.cca(X_FIT  = X_FIT[Fit_Index,],  Y_FIT  = Y_FIT[Fit_Index,],
                                           X_PRED = X_FIT[Pred_Index,], Y_PRED = Y_FIT[Pred_Index,], ncomp=ncomp,
                                           ProcrustX = CCA_OriginalData$xcoef, ProcrustY = CCA_OriginalData$ycoef)$variates)
      )

    }

    #Put cross-validated predictions back in original order and in a data frame
    Variates_CrossValidated = data.table::rbindlist(Variate_predictions)
    Variates_CrossValidated = Variates_CrossValidated[base::order(unlist(df_folds, FALSE, FALSE)),]

    #Bootstrap cross-validation predictions
    boot_i = base::sample(nrow(Variates_CrossValidated),nrow(Variates_CrossValidated), TRUE, NULL)
    Variates_CrossValidated_b = as.matrix(Variates_CrossValidated[boot_i,])

    #Estimate canonical correlations
    cc_CVboot[[b]] = Rfast::corpairs(as.matrix(Variates_CrossValidated_b[,1:ncomp]),as.matrix(Variates_CrossValidated_b[,(ncomp+1):(2*ncomp)]))
    # cc_CVboot = c(cc_CVboot , list(Rfast::corpairs(as.matrix(Variates_CrossValidated_b[,1:ncomp]),as.matrix(Variates_CrossValidated_b[,(ncomp+1):(2*ncomp)])))) #Not necessary

    #Estimate R2 for all outcome variables (with boot)
    Variates_CV_Scaled = as.matrix(Rfast::standardise(Variates_CrossValidated_b, center = TRUE, scale = TRUE))
    Y_FIT_Scaled =       as.matrix(Rfast::standardise(Y_FIT[boot_i,],center = TRUE, scale = TRUE))

    R2_matrix[[b]] =
      sapply(1:ncomp, function(ncomp_i)
        sapply(1:ncol(Y_FIT), function(y_i)
          R2quickcalc(X=Variates_CV_Scaled[,1:ncomp_i],Y=Y_FIT_Scaled[,y_i])
        )
      )

    #Standard Cross-Validation canonical correlation
    cc_CV[[b]] = Rfast::corpairs(as.matrix(Variates_CrossValidated[,1:ncomp]),as.matrix(Variates_CrossValidated[,(ncomp+1):(2*ncomp)]))


    if (UseProgressBar){
      utils::setTxtProgressBar(pb, value=b)
    }
  }


  # R2 Means
  R2_matrix = apply(base::simplify2array(R2_matrix), 1:2, mean)
  #Prettify output
  R2_matrix = data.frame(R2_matrix)
  rownames(R2_matrix) = colnames(Y_FIT)
  colnames(R2_matrix) = paste0("NVariates_",1:ncomp)

  # Cross-Validation Quantiles
  cc_CV_quantiles = do.call("rbind.data.frame", cc_CV)
  colnames(cc_CV_quantiles) = paste0("cc",1:ncomp)
  cc_CV_quantiles = apply(cc_CV_quantiles, 2, function(x) stats::quantile(x, probs = c(0.025,.5,.975), type=6))
  # cc_CV_pval = apply(cc_CV_quantiles, 2, function(x) boot_pval(x))

  # Cross-validation + Bootstrap Quantiles
  cc_CVBoot_quantiles = do.call("rbind.data.frame", cc_CVboot)
  colnames(cc_CVBoot_quantiles) = paste0("cc",1:ncomp)
  cc_CVBoot_quantiles2 = apply(cc_CVBoot_quantiles, 2, function(x) stats::quantile(x, probs = c(0.025,.5,.975), type=6))
  cc_CVBoot_pval = apply(cc_CVBoot_quantiles, 2, function(x) boot_pval(x))


  return(list(
    R2_matrix = R2_matrix,
    CrossValidationQuantiles = cc_CV_quantiles,
    CrossValidationBootstrapQuantiles = cc_CVBoot_quantiles2,
    CrossValidationBootstrapPvalues = cc_CVBoot_pval

  ))


}

