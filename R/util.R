EQM_cal <- function(i, bspline, data, allposknots, nknots, degree,
                    b.knots) {
    if ( isTRUE(bspline) ) 
              fit=lm(y ~ bs(x,knots=c(allposknots[2:(nknots+1),i]),
                                degree=degree), data = data)
          else 
              fit=lm(y ~ ns(x,knots=c(allposknots[2:(nknots+1),i]),
                                Boundary.knots=b.knots), data = data)
    mean(fit$residuals^2)
}
