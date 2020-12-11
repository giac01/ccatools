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
  res_ccatools = gb_CCA(X,Y,ncomp=NULL)


  #Check Canonical Correlations
  expect_equal(abs(res_base_r$cor),abs(res_ccatools$cc_fit))
  expect_equal(abs(res_base_r$cor),abs(res_ccatools$cc_pred))

  #Check Raw Coefficients
  expect_equal(abs(apply(res_base_r$xcoef,2, function(x) x/sqrt(sum(x^2)))),abs(res_ccatools$XLoadings))
  expect_equal(abs(apply(res_base_r$ycoef,2, function(x) x/sqrt(sum(x^2)))),abs(res_ccatools$YLoadings))

})

test_that("Check coef_boot ",{

  utils::data(iris)
  X = apply(iris[,1:2],2, as.numeric)
  Y = apply(iris[,3:4],2, as.numeric)

  res = gb_CCA(X,Y, ncomp=NULL)
  res_boot = coef_boot(X,Y,ncomp=2, Nboot=6000)

  #Check loadings look okay
  expect_equal(abs(res_boot$XLoadings_Quantiles[[1]]$original),abs(res$XLoadings[,1]))

  expect_equal(round(res_boot$XLoadings_Quantiles[[1]][,1]*100), c(90,-43) )
  expect_equal(round(res_boot$XLoadings_Quantiles[[1]][,3]*100), c(94,-35) )

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


