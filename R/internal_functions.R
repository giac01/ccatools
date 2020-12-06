
# --------------------------------------------------------------------------------------------------------------------------
# Internal Functions  ------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------

# Standard CCA function with additional functionality  ----------------------------------------------------------------------
# The function can fit a CCA model in one dataset, and use the loadings (coefficients) generated from this to predict CCA variate (latent variable) scores in a new dataset.
# In addition, the function can rotate the generated CCA loadings to a target matrix of X and Y loadings (ProcrustX and ProcrustY)
# Note: assumes NO missing data

gb_CCA = function(X_FIT,Y_FIT,X_PRED=NULL,Y_PRED=NULL,
                  ncomp=10,
                  ProcrustX = NULL, ProcrustY = NULL,
                  SafetyChecks=FALSE){
  # browser()
  #Check some basic things
  if (SafetyChecks){
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
  X = as.matrix(dataPreparation::fastScale(X_FIT, verbose=FALSE))
  Y = as.matrix(dataPreparation::fastScale(Y_FIT, verbose=FALSE))

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
  Yloadings = apply(as.matrix(Eig_V$vectors[,1:ncomp]),2,as.numeric)
  Xloadings = apply(as.matrix(Eig_W$vectors[,1:ncomp]),2,as.numeric)

  # Optional Rotation of Loadings

  if (!is.null(ProcrustX)){
    Xloadings = as.matrix(MCMCpack::procrustes(Xloadings  , ProcrustX)$X.new)

  }
  if (!is.null(ProcrustY)){
    Yloadings = as.matrix(MCMCpack::procrustes(Yloadings  , ProcrustY)$X.new)
  }

  # Estimate Canonical Variates

  Variates_Y =  Y_PRED %*% Yloadings
  colnames(Variates_Y) = paste0("Y", 1:ncomp)

  Variates_X =  X_PRED %*% Xloadings
  colnames(Variates_X) = paste0("X", 1:ncomp)

  Variates = cbind.data.frame(Variates_X, Variates_Y)

  cc_pred = Rfast::corpairs(Variates_X, Variates_Y)

  # Return Output

  out = list(
    #Loadings
    YLoadings = Yloadings,
    XLoadings = Xloadings,

    #Variates
    Variates = Variates,
    Variates_X = Variates_X,
    Variates_Y = Variates_Y,

    #Canonical Correlations
    cc_fit = cc,                            #Canonical correlations estimated from the input data (fitted canonical correlations)
    cc_pred = cc_pred                       #Canonical correlations estimated from the prediction data (predicted canonical correlations)
  )


  return(out)

}


# Function for VERY QUICK linear modelling  ---------------------------------------------------------------------------------------
# returns an R-squared as output
# WARNING - requires input to be scaled matrices!!
R2quickcalc = function(X,Y){
  # browser()
  out1= Rfast::lmfit(x=X,y=Y)
  Y = as.numeric(Y)
  return(1- sum(out1$residuals^2)/sum((Y-mean(Y))^2))
}

#Slower version, but does not require full rank matrices?
R2quickcalc_v2 = function(X,Y){
  # browser()
  out1= .lm.fit(x=X,y=Y)
  Y = as.numeric(Y)
  return(1- sum(out1$residuals^2)/sum((Y-mean(Y))^2))
}

# Estimate p=value from Bootstrap distribution ---------------------------------------------------------------------------------
boot_pval = function(x, null_val=0){
  x = na.omit(x)
  perc = length(which(x<null_val))/length(x)
  p_val = 1-abs(.50-perc)*2
  return(p_val)
}

# Function that is purely for simulation purposes ---------------------------------------------------------------------------------
# Automatically splis data in half - fits data in train half, and finds canonical corrlations in test half

gb_CCA_SplitHalfSim = function(X_FIT,Y_FIT, ncomp=10, ProcrustX = NULL, ProcrustY = NULL){

  df_folds = caret::createFolds(1:nrow(X_FIT), k=2)  #Split Data into halfs

  CCA_results = gb_CCA(X_FIT=X_FIT[df_folds[[1]],], Y_FIT=Y_FIT[df_folds[[1]],], X_PRED=X_FIT[df_folds[[2]],], Y_PRED=Y_FIT[df_folds[[2]],], ncomp=ncomp, ProcrustX = ProcrustX, ProcrustY = ProcrustY)

  CanonicalCorrelations = CCA_results$cc_pred # Estimated Canonical Correlations in hold-out dataset

  SampleSize = nrow(X_FIT[df_folds[[2]],])

  CC_ConfidenceIntervals = t(sapply(CanonicalCorrelations, function(x) CorrelationCIEstimator(x, n=SampleSize, alpha = .05)$CI))
  rownames(CC_ConfidenceIntervals) = paste0("cc",1:nrow(CC_ConfidenceIntervals))
  colnames(CC_ConfidenceIntervals) = c("CI_LB","CI_UB")
  CC_pvalues =  t(sapply(CanonicalCorrelations, function(x) CorrelationCIEstimator(x, n=SampleSize, alpha = .05)$p))

  return(list(
    CanonicalCorrelations=CanonicalCorrelations,
    CC_ConfidenceIntervals=CC_ConfidenceIntervals,
    CC_pvalues=CC_pvalues
  ))

}


# Pearson Correlation Confidence Interval Calculator  ---------------------------------------------------------------------------------

CorrelationCIEstimator = function(r, n, alpha=0.05){
  Fr = atanh(r)                 # Fisher Z Transform
  SE = 1/((n-3)^.5)             # Standard Error
  CI = qnorm(c(alpha/2,1-alpha/2), mean=Fr, sd=SE)
  CI = tanh(CI)
  # p  = (1-pnorm(abs(Fr), mean=0, sd=SE))*2    # Fisher Z P value
  t = r*sqrt((n-2)/(1-r^2))       # P-value estimated from t-distribution
  p  = (1-pt(abs(t), df=n-2))*2
  return(list(CI=CI,p=p))
}
