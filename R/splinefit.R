library(lubridate)
library(dplyr)
library(ggplot2)
library(splines)
library(stringr)
library(parallel)
library(Rcpp)

#' @title Automatic spline-based model selection
#' @description Provides optimal number of internal knots and their locations in spline-based models.
#' @param x Vector of function arguments.
#' @param y Vector of noisy function values.
#' @param nlimknot Maximum number of internal knots.
#' @param lambda Tunning parameter of the penalized least squares estimation process.
#' @param delta Minimal distance value between two consecutive knots.
#' @param degree Degree of the spline basis.
#' @param bspline External knots vector. Default is the (min(x),max(x)).
#' @param b.knots Logical value. If TRUE, B-spline is fitted. Else, natural splines.
#'
#' @return
#' @export
#'
#' @examples
 splinefit = function(x,y,nlimknot,lambda,delta,degree,bspline=TRUE,
                               b.knots=NULL){
  
  fixknots = matrix(0,nlimknot+1,nlimknot+1)
  if ( is.null(b.knots) ) {
    if ( isTRUE(bspline) ) {
      first.knot = sort(x)[round(length(x)*0.25)]
      last.knot = max(x)
      b.knots = c(first.knot, last.knot)
      print(b.knots)
    }
    else {
      first.knot = sort(x)[round(length(x) *
                                               0.10)] 
      last.knot = sort(x)[round(length(x) *
                                              0.95)]
      b.knots = c(first.knot, last.knot)
    }
  }
  
  EQMpen = vector()  
  for(i in 1:(nlimknot+1)){
    a = knotfixfit(x,y,i-1,delta,degree,b.knots,bspline)
    EQMpen[i] = a$min.eqm + lambda*i
    fixknots[i,1:length(a$knots)] = a$knots
  }
  
  min.nknot = which.min(EQMpen)
  if ( min.nknot == 1 ) {bestknots = NULL
  }else {bestknots = fixknots[which.min(EQMpen),1:(min.nknot-1)]}
  
  if ( isTRUE(bspline) ){
    bestfitall=lm(y ~ bs(x,knots = c(bestknots),degree=degree))
  }else{
    bestfitall=lm(y ~ ns(x,knots = c(bestknots),Boundary.knots=b.knots))}
  
  return(list(bestknots, bestfitall))
}
