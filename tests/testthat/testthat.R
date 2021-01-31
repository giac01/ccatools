library(testthat)
library(ccatools)

# devtools::load_all()

test_that("Internal Package CCA function gives identical output as stats::cancor ",{

  set.seed(100)
  X = scale(sapply(1:50, function(x) stats::rnorm(200)))
  Y = scale(sapply(1:50, function(x) stats::rnorm(200)))

  #Get base R result
  res_base_r = stats::cancor(X,Y)

  #Get package R result
  res_ccatools = .cca(X,Y,ncomp=NULL)


  #Check Canonical Correlations
  expect_equal(abs(res_base_r$cor),abs(res_ccatools$cc_fit))
  expect_equal(abs(res_base_r$cor),abs(res_ccatools$cc_pred))

  #Check Raw Coefficients
  expect_equal(abs(apply(res_base_r$xcoef,2, function(x) x/sqrt(sum(x^2)))),abs(res_ccatools$xcoef))
  expect_equal(abs(apply(res_base_r$ycoef,2, function(x) x/sqrt(sum(x^2)))),abs(res_ccatools$ycoef))

})


test_that("Test functions work when there is only one outcome",{

  set.seed(100)
  X = scale(sapply(1:50, function(x) stats::rnorm(200)))
  Y = as.matrix(scale(sapply(1:1, function(x) stats::rnorm(200))))

  #Get base R result

  #Get package R result
  res_ccatools = .cca(X,Y,ncomp=NULL)
  expect_error(
    cca_splithalf(X_FIT = X[1:100,],    Y_FIT  = as.matrix(Y[1:100,]),
                  X_PRED = X[101:200,], Y_PRED = as.matrix(Y[101:200,]),
                  alpha = .05
    ), NA
  )

})

test_that("Check coef_boot ",{

  utils::data(iris)
  X = apply(iris[,1:2],2, as.numeric)
  Y = apply(iris[,3:4],2, as.numeric)

  res = .cca(X,Y, ncomp=NULL)
  res_boot = coef_boot(X,Y,ncomp=2, Nboot=10000)

  #Check loadings look okay
  expect_equal(abs(res_boot$xcoef_Quantiles[[1]]$original),abs(res$xcoef[,1]))

  expect_equal(round(res_boot$xcoef_Quantiles[[1]][,1]*100), c(90,-43) )
  expect_equal(round(res_boot$xcoef_Quantiles[[1]][,3]*100), c(94,-35) )

})

test_that("Check cca_splithalf ",{

  utils::data(iris)
  X = apply(iris[,1:2],2, as.numeric)
  Y = apply(iris[,3:4],2, as.numeric)

  res = cca_splithalf(X_FIT = X[1:75,],    Y_FIT  = Y[1:75,],
                      X_PRED = X[76:150,], Y_PRED = Y[76:150,],
                      alpha = .05
                      )

  #Check Results look okay
    expect_equal(res$predicted_cc, c(0.798424615016857, 0.131374767458812))

    expect_equal(c(t(res$confint_cc)), c(0.697937337494885, -0.0985247974926762,0.868080011288514,0.347961426057066 ))

    expect_equal(c(t(res$pvalue_cc)), c(0, 0.261222183241876))

})

test_that("Check R2 Estimates Are Good ",{

  N = 1000
  x = rnorm(N, mean= 10)

  x_mat = sapply(1:10, function(i) x + 1*rnorm(N))
  y = x + rnorm(N)
  dat = cbind.data.frame(x_mat,y)
  out = summary(lm(y ~ ., data=data.frame(scale(dat))))


  expect_equal(
    .R2quickcalc_v2(scale(x_mat),scale(y)),
    .R2quickcalc(scale(x_mat),scale(y)),
    out$r.squared


  )
})

## ---- Ignore below ----

## Code to test function
#
# rm(list=ls())
# set.seed(100)
# X = scale(sapply(1:30, function(x) stats::rnorm(1000)))
# Y = scale(sapply(1:22, function(x) stats::rnorm(1000)))
#
# X_FIT = X[1:500,]
# Y_FIT  = as.matrix(Y[1:500,])
# X_PRED = X[501:1000,]
# Y_PRED = as.matrix(Y[501:1000,])
# alpha = .05
# ProcrustX = NULL
# ProcrustY = NULL
# ncomp=10
# rm(X,Y)

# #
# rm(list=ls())
# set.seed(100)
# N = 200 # MUST BE EVEN NUMBER
# x = rnorm(N)
# y = scale(rnorm(N) + x*0)
# X = scale(sapply(1:30, function(z) x + stats::rnorm(N)))
# Y = scale(sapply(1:22, function(z) y + stats::rnorm(N)))
#
# X_FIT = X[1:(N/2),]
# Y_FIT  = as.matrix(Y[1:(N/2),])
# X_PRED = X[(N/2+1):N,]
# Y_PRED = as.matrix(Y[(N/2+1):N,])
# alpha = .05
# ProcrustX = NULL
# ProcrustY = NULL
# ncomp=10
# rm(X,Y,x,y)
#
# model_test = ccatools::cca_splithalf(X_FIT = X_FIT,
#                                Y_FIT = Y_FIT,
#                                X_PRED = X_PRED,
#                                Y_PRED = Y_PRED,
#                                alpha = .05,
#                                ProcrustX = NULL,
#                                ProcrustY = NULL,
#                                ncomp=1
# )
#
#
#
# # model_test$R2_matrix_unbiased[,-11]>model_test$R2_matrix_unbiased[,11]
# model_test$R2_matrix_unbiased
#
# model_test$R2_matrix

