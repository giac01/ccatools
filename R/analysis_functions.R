#' Split-Half CCA code
#'
#' Run CCA model in training dataset, and and validate performance in testing dataset.
#' The function estimates confidence intervals and p-values using the ... .
#'
#' @param X_FIT Numeric Matrix [N, P1] containing the training dataset predictor variables.
#' @param Y_FIT Numeric Matrix [N, P2] containing the training dataset outcome variables.
#' @param X_PRED Numeric Matrix [N, P1] containing the testing dataset predictor variables. Variables should be ordered in the same way as for X_FIT.
#' @param Y_PRED Numeric Matrix [N, P1] containing the testing dataset outcome variables. Variables should be ordered in the same way as for Y_FIT.
#' @param ncomp Numeric Scalar. Number of CCA components to keep in analyses. Must be equal to or less than min(P1,P2).
#' @param alpha Numeric Scalar. Alpha level for estimating a 100(1-alpha)\% confidence interval for each canonical correlation. Default is .05 for estimating a 95\% confidence interval.
#'
#' @export
#'
cca_splithalf = function(X_FIT,Y_FIT,X_PRED,Y_PRED,
                            ncomp=NULL, alpha = 0.05){

  model_results = gb_CCA(X_FIT=X_FIT,Y_FIT=Y_FIT,X_PRED=X_PRED,Y_PRED=Y_PRED,
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
#' This function runs the gb_CCA Canonical Correlation Analysis function multiple times to assess variability in the CCA loadings and canonical correlations
#' Because bootstrap resampling can change the order of canonical variates that are extracted, or sign flipping can occur
#' in some cases (i.e. a very similar latent variable is extracted but on some occasions the loadings are mostly positive or negative), we rotate the loadings
#' in each bootstrap resample to map onto the loadings generated from the full, raw input datsets.
#'
#'
#' @param X_FIT Numeric Matrix [N, P1] containing the predictor variables.
#' @param Y_FIT Numeric Matrix [N, P2] containing the outcome variables.
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

  #Run model on full dataset
  CCA_OriginalData = gb_CCA(X_FIT=X_FIT, Y_FIT=Y_FIT, X_PRED=NULL, Y_PRED=NULL,
                            ProcrustX = ProcrustX, ProcrustY = ProcrustY,
                            ncomp=ncomp)

  #Lists to store bootstrap data in
  YLoadings_ROTATED = list()  #Store rotated loadings from bootstrap resamples in here
  XLoadings_ROTATED = list()
  cc = list()                #Store canonical correlations from bootstrap resamples in hre

  for(i in 1:Nboot){
    # Randomly bootstrap resample data and fit CCA model.
    BootResample = base::sample(nrow(X_FIT), replace = TRUE)


    #Run CCA on bootstrap resampled data, rotating loadings using procrustes to map  onto the original model
    CCA_BootData = gb_CCA(X_FIT=X_FIT[BootResample,], Y_FIT=Y_FIT[BootResample,],
                          ProcrustX = CCA_OriginalData$XLoadings, ProcrustY = CCA_OriginalData$YLoadings,
                          X_PRED=NULL, Y_PRED=NULL,
                          ncomp=ncomp)

    XLoadings_ROTATED[[i]] = CCA_BootData$XLoadings
    YLoadings_ROTATED[[i]] = CCA_BootData$YLoadings
    cc[[i]]  = CCA_BootData$cc_pred

    utils::setTxtProgressBar(pb, value=i)

  }
  #Estimate quantiles from boostrap distribution

  XLoadings_Quantiles = apply(base::simplify2array(XLoadings_ROTATED), 1:2, stats::quantile, prob = c(0.025, .5, 0.975))
  XLoadings_Quantiles = lapply(1:ncomp, function(i) XLoadings_Quantiles[,,i])
  for(i in 1:ncomp){
    colnames(XLoadings_Quantiles[[i]]) = colnames(X_FIT)
    XLoadings_Quantiles[[i]] = data.frame(t(XLoadings_Quantiles[[i]]))
    XLoadings_Quantiles[[i]]$original = CCA_OriginalData$XLoadings[,i]
  }

  YLoadings_Quantiles = apply(base::simplify2array(YLoadings_ROTATED), 1:2, stats::quantile, prob = c(0.025, .5, 0.975))
  YLoadings_Quantiles = lapply(1:ncomp, function(i) YLoadings_Quantiles[,,i])
  for(i in 1:ncomp){
    colnames(YLoadings_Quantiles[[i]]) = colnames(Y_FIT)
    YLoadings_Quantiles[[i]] = data.frame(t(YLoadings_Quantiles[[i]]))
    YLoadings_Quantiles[[i]]$original = CCA_OriginalData$YLoadings[,i]
  }

  cc_quantiles = apply(base::simplify2array(cc), 1, stats::quantile, prob = c(0.025, .5, 0.975))


  # Ouput Data

  out = list(
    XLoadings_Quantiles=XLoadings_Quantiles,
    YLoadings_Quantiles=YLoadings_Quantiles,
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
#' @param X_FIT Numeric Matrix [N, P1] containing the predictor variables.
#' @param Y_FIT Numeric Matrix [N, P2] containing the outcome variables.
#' @param ncomp Numeric Scalar. Number of CCA components to keep in analyses. Must be equal to or less than min(P1,P2).
#' @param Nboot Numeric Scalar. Number of times to repeat k-fold cross-validation (yes its confusing it says "boot").
#' @param Nfolds Numeric Scalar. Number of
#' @param UseProgressBar Logical. Whether to show progress bar.
#' @param UseProcrustes  Logical. Whether to use procrustes analysis.
#'
#' @export
#'
cca_cv_boot = function(X_FIT,Y_FIT, ncomp=10, Nboot=30, Nfolds=10,ProcrustX = NULL, ProcrustY = NULL, UseProgressBar=TRUE, UseProcrustes=TRUE){
  # browser()
  if (UseProgressBar){
    pb <- utils::txtProgressBar(min = 0, max = Nboot, style = 3)
  }

  CCA_OriginalData = gb_CCA(X_FIT=X_FIT, Y_FIT=Y_FIT, X_PRED=NULL, Y_PRED=NULL, ncomp=ncomp, ProcrustX = ProcrustX, ProcrustY = ProcrustY)

  if (UseProcrustes==FALSE){
    CCA_OriginalData$XLoadings = NULL # By setting this to NULL, the gb_CCA function called below will not use procrustes rotations
    CCA_OriginalData$YLoadings = NULL
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
                               list(gb_CCA(X_FIT  = X_FIT[Fit_Index,],  Y_FIT  = Y_FIT[Fit_Index,],
                                           X_PRED = X_FIT[Pred_Index,], Y_PRED = Y_FIT[Pred_Index,], ncomp=ncomp,
                                           ProcrustX = CCA_OriginalData$XLoadings, ProcrustY = CCA_OriginalData$YLoadings)$Variates)
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

