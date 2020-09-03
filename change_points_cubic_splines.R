####################################################################################
###### Obtem melhor ajuste para os dados do covid para um modelo de
###### regressão com splines cúbicos naturais, com número e posições
###### variáveis para os knots, utilizando métodos de seleção de
###### modelos e detecção de pontos de mudança
######
###### Autores: Alex Rodrigo, Florencia Leonardi e Magno Severino
###### Instituto de Matemática e Estatística - Universidade de São Paulo
###### Data: 03/06/2020
####################################################################################


library(lubridate)
library(dplyr)
library(ggplot2)
library(splines)
library(stringr)
library(parallel)
library(Rcpp)


sourceCpp("not_able_knots.cpp")

##### INÍCIO FUNÇÕES

##### Best fit for a fixed number of knots  ######

splineSelection.fitfix = function(covdata,nknots,delta,degree,
                           b.knots=c(min(covdata$index),max(covdata$index)),
                           bspline=TRUE){
  
  if( (b.knots[2] - b.knots[1]) / delta < (nknots+1) ){
    print("nknots is too large with respect to delta: returning NA for this number of nknots\n")
    return(list(knots=NA,min.eqm=NA,model=NULL))
  }
  
  posknots <- covdata$index[!(covdata$index %in% c(min(covdata$index), max(covdata$index)))] %>% 
    as.integer() %>% sort()

  posknots <- posknots[posknots >= b.knots[1] & posknots <= b.knots[2]]
  
  if( nknots > 0){

      allposknots = rbind(min(covdata$index),
                         combn(posknots,nknots),
                         max(covdata$index))


      notAble <- not_able_knots(allposknots, delta)
      allposknots = allposknots[, -notAble]
      ## allposknots = allposknots[,apply(apply(allposknots,FUN=dist,2),
      ## FUN=min,2)>delta]

      if( ncol(allposknots) == 0){
          print("allposknots with 0 columns: decrease nknots or delta, or aumentgs sample size\n")
          return(list(knots=NA,min.eqm=NA,model=NULL))
      }
      
      EQM = NA

      ret = mclapply(1:ncol(allposknots), EQM_cal, bspline, covdata,
                     allposknots, nknots, degree, b.knots,
                     mc.cores=detectCores())
      EQM = simplify2array(ret)      

      ## for(i in 1:ncol(allposknots)){
      ##     if ( isTRUE(bspline) ) 
      ##         fit=lm(cases ~ bs(index,knots=c(allposknots[2:(nknots+1),i]),
      ##                           degree=degree), data = covdata )
      ##     else 
      ##         fit=lm(cases ~ ns(index,knots=c(allposknots[2:(nknots+1),i]),
      ##                           Boundary.knots=b.knots), data = covdata )
      ##     EQM[i] = mean(fit$residuals^2)
      ## }

      chosenknots = allposknots[2:(nknots+1),which.min(EQM)]
      
      if ( isTRUE(bspline) )
          bestfit=lm(cases ~ bs(index,knots = c(chosenknots),degree=degree),
                     data = covdata )
      else
          bestfit=lm(cases ~ ns(index,knots = c(chosenknots),
                                Boundary.knots=b.knots), data = covdata ) 
  }
  else {
      chosenknots = NA
      if ( isTRUE(bspline) )
          bestfit=lm(cases ~ bs(index,knots=NULL,degree=degree),
                     data = covdata )
      else
          bestfit=lm(cases ~ ns(index,knots=NULL,Boundary.knots=b.knots),
                     data = covdata ) 
      EQM = mean(bestfit$residuals^2)
  } 
    
    return(list(knots=chosenknots,min.eqm=min(EQM),model=bestfit))
    
}

#### Best fit for all number of knots (until a limit nlimknot) ####

splineSelection.fit = function(covdata,nlimknot,lambda,delta,degree,bspline=TRUE,
                        b.knots=NULL,nome='nome'){
    
    fixknots = matrix(0,nlimknot+1,nlimknot+1)
    if ( is.null(b.knots) ) {
        if ( isTRUE(bspline) ) {
            first.knot = sort(covdata$index)[round(length(covdata$index) *
                                                   0.25)]
            last.knot = max(covdata$index)
            b.knots = c(first.knot, last.knot)
            print(b.knots)
        }
        else {
            first.knot = sort(covdata$index)[round(length(covdata$index) *
                                                   0.10)] 
            last.knot = sort(covdata$index)[round(length(covdata$index) *
                                                  0.95)]
            b.knots = c(first.knot, last.knot)
        }
    }
    
    EQMpen = vector()  
    for(i in 1:(nlimknot+1)){
        a = splineSelection.fitfix(covdata,i-1,delta,degree,b.knots,bspline)
        EQMpen[i] = a$min.eqm + lambda*i
        fixknots[i,1:length(a$knots)] = a$knots
    }
    
    min.nknot = which.min(EQMpen)
    if ( min.nknot == 1 ) {bestknots = NULL
    }else {bestknots = fixknots[which.min(EQMpen),1:(min.nknot-1)]}
    
    if ( isTRUE(bspline) ){
        bestfitall=lm(cases ~ bs(index,knots = c(bestknots),degree=degree), data = covdata )
    }else{
        bestfitall=lm(cases ~ ns(index,knots = c(bestknots),Boundary.knots=b.knots), data = covdata ) }
    
    return(list(bestknots, bestfitall))
}

## Calcula EQM da i-ésima combinação 
EQM_cal <- function(i, bspline, covdata, allposknots, nknots, degree,
                    b.knots) {
    if ( isTRUE(bspline) ) 
              fit=lm(cases ~ bs(index,knots=c(allposknots[2:(nknots+1),i]),
                                degree=degree), data = covdata )
          else 
              fit=lm(cases ~ ns(index,knots=c(allposknots[2:(nknots+1),i]),
                                Boundary.knots=b.knots), data = covdata )
    mean(fit$residuals^2)
}
