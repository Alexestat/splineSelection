knotfixfit = function(x,y,nknots,delta,degree,
                                  b.knots=c(min(x),max(x)),
                                  bspline=TRUE){
  
  if( (b.knots[2] - b.knots[1]) / delta < (nknots+1) ){
    print("nknots is too large with respect to delta: returning NA for this number of nknots\n")
    return(list(knots=NA,min.eqm=NA,model=NULL))
  }
  
  posknots <- x[!(x %in% c(min(x), max(x)))] %>% 
    as.integer() %>% sort()
  
  posknots <- posknots[posknots >= b.knots[1] & posknots <= b.knots[2]]
  
  if( nknots > 0){
    
    allposknots = rbind(min(x),
                        combn(posknots,nknots),
                        max(x))
    
    
    notAble <- not_able_knots(allposknots, delta)
    allposknots = allposknots[, -notAble]
    
    if( ncol(allposknots) == 0){
      print("allposknots with 0 columns: decrease nknots or delta, or aumentgs sample size\n")
      return(list(knots=NA,min.eqm=NA,model=NULL))
    }
    
    EQM = NA
    
    ret = mclapply(1:ncol(allposknots), EQM_cal, bspline,data.frame(x,y),
                   allposknots, nknots, degree, b.knots,
                   mc.cores=detectCores())
    EQM = simplify2array(ret)      
    
    chosenknots = allposknots[2:(nknots+1),which.min(EQM)]
    
    if ( isTRUE(bspline) )
      bestfit=lm(y ~ bs(x,knots = c(chosenknots),degree=degree))
    else
      bestfit=lm(y ~ ns(x,knots = c(chosenknots),
                        Boundary.knots=b.knots)) 
  }
  else {
    chosenknots = NA
    if ( isTRUE(bspline) )
      bestfit=lm(y ~ bs(x,knots=NULL,degree=degree))
    else
      bestfit=lm(y ~ ns(x,knots=NULL,Boundary.knots=b.knots)) 
    EQM = mean(bestfit$residuals^2)
  } 
  
  return(list(knots=chosenknots,min.eqm=min(EQM),model=bestfit))
  
}