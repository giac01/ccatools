.sortVar = function(dat){

  cormat = 1 - base::abs(stats::cor(dat, use="pairwise.complete.obs"))
  cormat = stats::as.dist(cormat)
  clustres = stats::hclust(cormat)

  return(clustres$order)
}




# Confidence Interval Function Required for PlotCorrelationMatrix
.R_ConInt = function(x,y){
  # browser()
  if(base::is.na(stats::cor(x,y, use="pairwise.complete.obs"))){
    # return(list(CI="0", N=nrow(na.omit(cbind.data.frame(x,y))), CI_l="0"))
    return(base::list(CI="0", N="0", CI_l="0"))

  } else {
    CI = stats::cor.test(x,y)$conf.int[1:2]
    CI = base::gsub("^\\s","",gsub("0\\.", "\\.",base::format(CI, digits=0, nsmall=2, perl=FALSE)))
    N = base::as.numeric(cor.test(x,y)$parameter + 2)
    CI_l =     base::paste0( "[",CI[1],", ",CI[2],"]", sep="")
    out =    base::paste0(N,"\n[",CI[1],", ",CI[2],"]", sep="")
    return(list(CI=out, N=N, CI_l=CI_l))
  }
}



#' Create a ggplot2 Correlation Matrix
#'
#' @Decription
#'
#' The function provides a wrapper around ggplot2 to quickly create a correlation matrix.
#'
#' @param dat Input dataframe or matrix (*do not input a correlation matrix*).
#' @param Variables_Labels Character vector of variable Labels, corresponding to each column in dat. If missing (NULL) then colnames(dat) will be used.
#' @param textadjust Scalar. Adjust text size by a magnification factor.
#' @param includeN Logical. Include sample size on upper diagononal (TRUE) or leave blank (FALSE).
#' @param reportCI Logical. Include confidence interval on upper diagonal (TRUE) or leave blank (FALSE).
#' @param low_colour Logical. Hex colour code for the lowest correlation.
#' @param high_colour Logical. Hex colour code for the highest correlation.
#' @param abs_colour Logical. If TRUE, will use the absolute correlation (i.e. ignoring whether the correlation is positive or negative) for determining square colour.
#'
#' @return ggplot
#'
#' @examples
#' X = sapply(1:10, function(i) rnorm(100))
#' X = as.data.frame(X)
#' My_Labels = c(paste0("Predictor ",1:5), paste0("Outcome ",1:5))
#'
#' plotcor(X, Variables_Labels = My_Labels, includeN = TRUE, reportCI = FALSE)
#'
#' @export
#'
plotcor = function(dat,
                                 Variables_Labels=NULL, textadjust=2, includeN=TRUE, reportCI=TRUE,
                                 low_colour="#EDCB64", high_colour="#B62A3D", abs_colour=TRUE,
                                 cluster_variables = FALSE
                                 ){
  if (!base::is.data.frame(dat)) {dat=base::as.data.frame(dat)}

  if (cluster_variables) {
    new_order = ccatools:::.sortVar(dat)
    dat = dat[,new_order]
    if (!base::is.null(Variables_Labels)){
      Variables_Labels = Variables_Labels[new_order]
    }
  }

  Variables = base::colnames(dat)
  if(is.null(Variables_Labels)){
    Variables_Labels = base::colnames(dat)
  }

  matrix_scores = dat
  Mat_Cor =gsub("^0","",gsub("^ +","",gsub("^-0","-", format(cor(matrix_scores, use="pairwise.complete.obs"), digits=0, nsmall=2)))) #Correlation matrix
  if(abs_colour){
    Mat_Cor_fill = base::abs(stats::cor(matrix_scores, use="pairwise.complete.obs"))  #Correlation matrix for table fill
  } else {
    Mat_Cor_fill = (stats::cor(matrix_scores, use="pairwise.complete.obs"))  #Correlation matrix for table fill
  }

  Mat_Cor_fill[base::lower.tri(Mat_Cor_fill,diag = TRUE)]=NA

  #Matrix on Ns per comparison - lower triag
  Mat_N = sapply(Variables, function(x)
    sapply(Variables, function(y)
      base::nrow(stats::na.omit(base::data.frame(dat[,base::unique(c(x,y))])))
    ))

  #Confidence interval information
  if(reportCI){
    Mat_CI = sapply(Variables, function(x)
      sapply(Variables, function(y)
        ccatools:::.R_ConInt(as.numeric(dat[,x]),base::as.numeric(dat[,y]))$CI
      ))
    base::diag(Mat_CI) = base::diag(Mat_N)
    Mat_N = Mat_CI
  }

  #Create Dataframe For ggplot to Read
  PlotMat = Mat_Cor
  if(includeN){
    PlotMat[base::lower.tri(PlotMat, diag=TRUE)]=Mat_N[base::lower.tri(Mat_N, diag=TRUE)]
  }
  if(!includeN){
    PlotMat[base::lower.tri(PlotMat, diag=TRUE)]=""
  }
  base::colnames(PlotMat) = Variables_Labels ;  base::rownames(PlotMat) = Variables_Labels

  PlotMat = base::data.frame(reshape2::melt(PlotMat), stringsAsFactors = FALSE)
  head(PlotMat)

  PlotMat$value = (base::as.character(PlotMat$value))
  PlotMat$ValueFill = base::as.numeric(t(c(Mat_Cor_fill)))
  PlotMat$Var2 = base::factor(PlotMat$Var2, levels=rev(base::levels(PlotMat$Var2)))


  OutPlot =
    ggplot2::ggplot(data = PlotMat, ggplot2::aes(x=Var1, y=Var2,fill=ValueFill))+# + geom_point(aes(size=value^2,alpha=value^4))+
    ggplot2::geom_tile() + ggplot2::labs(x=NULL, y=NULL) +
    ggplot2::theme(axis.text = ggplot2::element_text(size=5*textadjust)) +
    ggplot2::geom_text(ggplot2::aes(label=value), size=1.4*textadjust) +
    jtools::theme_apa() +
    #scale_fill_brewer(palette=1,na.value="grey")+
    ggplot2::scale_fill_continuous(na.value="gray94",low=low_colour, high=high_colour) +
    ggplot2::theme(legend.position = "#DAECED") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,hjust=0.95,vjust=0.2)) +
    ggplot2::coord_fixed()

  OutPlot

  return(OutPlot)
}

#' mapnornmal - Map any continuous variable onto a Normal Distribution
#'
#' Transforms any numeric vector onto a continuous distribution.
#'
#' The function preserves the original order of values, but shifts the location of values so that the overall distribution maps perfectly onto a normal distribution.
#'
#' This works by first finding the percentile ranks (i.e. rank individuals between 0 and 1), and then use the inverse cumulative distribution function for a standard normal distribution to convert those ranks to the equivalent position on a normal distribution.
#'
#' Two additional tricks are required.
#'
#' First, after ranking values from 1 to N (N = the number of values), we minus 1/2 then divide by N to get percentile ranks.
#' This is to avoid having percentile ranks of 0% or 100% - as they map onto -/+ infinity!
#'
#' Second problem is dealing with tied values. The approach here is to initially sort tied values randomly, assign them all a unique normal score, then average across the normal scores for each tied value.
#' This yields a different result than first averaging the percentile rank, and finding a normal score from that.
#'
#' After transforming percentile ranks to the normal distribution, additional scaling and centering is performed to make sure the vector mean is 0 with unit variance.
#'
#' Warning 1: The default behavior is that the function will not transform dichotomous variables. It's generally only advisable to use this function for continuous variables with minimal repeated values.
#'
#' Warning 2: For linear regression your outcome variable does not need to be normally distributed! It may be helpful to use this function for predictor variables, but transforming outcomes may be inadvisable....
#'
#' @param x Numeric Vector that your want to map onto a Normal distribution.
#' @param MinDim Numeric Scalar (default = 2). If the number of unique values in a vector is equal to or less than MinDim, then no the function will not transform the variable. This is because it is inadvisable to use this function on dichotmous variables.
#' @param rescale Logical (default = TRUE). If true, will scale() the output to ensure a 0 mean and unit variance.
#'
#' @return
#' @export
#'
#' @examples
map_normal = function(x, MinDim=2, rescale=TRUE){
    if (!is.numeric(x)) stop("Uh oh, make sure your input is a numeric vector!")

  #Ranks the kids on unique integers (for ties, randomly order these),then finds a average normal score
    percentilerank = (base::rank(x, na.last="keep", ties.method = "random")-.5)/(base::length(stats::na.omit(x)))
    normalscore = stats::qnorm(percentilerank)
  #average over z-scores with the same value.
    for(i in base::unique(x)){
      if(base::length(base::which(x==i))>1){
        normalscore[base::which(x==i)]=mean(normalscore[x==i],na.rm = TRUE)
      }
    }
    out = base::as.numeric(normalscore)

    if (rescale==TRUE){
      out = base::as.numeric(base::scale(out,center = TRUE, scale=TRUE))
    }

  #Minimum number of unique values in data allowed. This is so that we can avoid transforming binary data types or likert scales if we don't want them to be transformed
  ObservedNumberUniqueValues = base::length(base::unique(stats::na.omit(x)))
  if (ObservedNumberUniqueValues <= MinDim){
    out=x
  }

  return(out)
}


# Function required for concert
.swapsies = function(x, input, replacement){
  return(replacement[match(x, input)])
}

recode_check = function(x, input, output, verbose = TRUE){

  before=c(t(x))
  after = .swapsies(before, input, output)

  if (verbose){
    print(table(before))
    print(table(after))
  }

  return(after)

}




