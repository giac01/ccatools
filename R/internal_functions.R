.isBinaryVariable = function(x){

  return(2==base::length(base::table(x,useNA="no")))

}

# Converts into binary variable (0, 1)

.makeBinary = function(x){

  return(base::as.numeric(x>base::min(x, na.rm = TRUE)))

}

.inverseLogit = function(x){
  return(1/(1+exp(-1*(x))))
}


# Coefficient of Determination Estimation
## https://en.wikipedia.org/wiki/Coefficient_of_determination

.R2 = function(X,Y){
  if (nrow(X)!=nrow(Y)) stop(".R2 Fail; nrow don't match!!")
  Res = Y - X #Matrix of residuals
  Yt = base::apply(base::as.matrix(Y), 2, function(x) x - mean(x, na.rm=TRUE) )

  SSres =  base::crossprod(Res, Res)        # Residual Sums of Squares
  SStot =  base::crossprod(Yt, Yt)          # Y (total) Sums of Squares

  R2 = base::diag(1 - (SSres/SStot))

  return(R2)
}


# Add intercept to matrix
.addintercept = function(matrix){
  intercept = rep(1, nrow(matrix)) # add intercept column to training data

  return(
    (
      cbind(
        intercept,
        matrix
      )
    )
  )

}


# Function for VERY QUICK linear modelling  ---------------------------------------------------------------------------------------
# returns an R-squared as output
# WARNING - requires input to be scaled matrices!!
.R2quickcalc = function(X,Y){
  # browser()
  intercept = base::rep(1,nrow(X))
  X = base::cbind(intercept, X)
  out1= Rfast::lmfit(x=X,y=Y)
  Y = as.numeric(Y)
  return(1- sum(out1$residuals^2)/sum((Y-mean(Y))^2))
}

#Slower version, but does not require full rank matrices?
.R2quickcalc_v2 = function(X,Y){
  # browser()
  intercept = base::rep(1,nrow(X))
  X = base::cbind(intercept, X)
  out1= stats::.lm.fit(x=X,y=Y)
  Y = as.numeric(Y)
  return(1- sum(out1$residuals^2)/sum((Y-mean(Y))^2))
}

# Estimate p=value from Bootstrap distribution ---------------------------------------------------------------------------------
.boot_pval = function(x, null_val=0){
  x = stats::na.omit(x)
  perc = length(which(x<null_val))/length(x)
  p_val = 1-abs(.50-perc)*2
  return(p_val)
}

# Function that is purely for simulation purposes ---------------------------------------------------------------------------------
# Automatically splis data in half - fits data in train half, and finds canonical corrlations in test half

gb_CCA_SplitHalfSim = function(X_FIT,Y_FIT, ncomp=10, ProcrustX = NULL, ProcrustY = NULL){

  df_folds = caret::createFolds(1:nrow(X_FIT), k=2)  #Split Data into halfs

  CCA_results = .cca(X_FIT=X_FIT[df_folds[[1]],], Y_FIT=Y_FIT[df_folds[[1]],], X_PRED=X_FIT[df_folds[[2]],], Y_PRED=Y_FIT[df_folds[[2]],], ncomp=ncomp, ProcrustX = ProcrustX, ProcrustY = ProcrustY)

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

.CorrelationCIEstimator = function(r, n, alpha=0.05){
  Fr = base::atanh(r)                 # Fisher Z Transform
  SE = 1/((n-3)^.5)             # Standard Error
  CI = stats::qnorm(c(alpha/2,1-alpha/2), mean=Fr, sd=SE)
  CI = base::tanh(CI)
  # p  = (1-pnorm(abs(Fr), mean=0, sd=SE))*2    # Fisher Z P value
  t = r*base::sqrt((n-2)/(1-r^2))       # P-value estimated from t-distribution
  p  = (1-stats::pt(base::abs(t), df=n-2))*2
  return(base::list(CI=CI,p=p))
}
