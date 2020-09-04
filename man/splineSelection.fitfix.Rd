\name{splineSelection.fitfix}
\alias{splineSelection.fitfix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
splineSelection.fitfix(covdata, nknots, delta, degree, b.knots = c(min(covdata$index), max(covdata$index)), bspline = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covdata}{
%%     ~~Describe \code{covdata} here~~
}
  \item{nknots}{
%%     ~~Describe \code{nknots} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
}
  \item{b.knots}{
%%     ~~Describe \code{b.knots} here~~
}
  \item{bspline}{
%%     ~~Describe \code{bspline} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (covdata, nknots, delta, degree, b.knots = c(min(covdata$index), 
    max(covdata$index)), bspline = TRUE) 
{
    if ((b.knots[2] - b.knots[1])/delta < (nknots + 1)) {
        print("nknots is too large with respect to delta: returning NA for this number of nknots\n")
        return(list(knots = NA, min.eqm = NA, model = NULL))
    }
    posknots <- covdata$index[!(covdata$index \%in\% c(min(covdata$index), 
        max(covdata$index)))] \%>\% as.integer() \%>\% sort()
    posknots <- posknots[posknots >= b.knots[1] & posknots <= 
        b.knots[2]]
    if (nknots > 0) {
        allposknots = rbind(min(covdata$index), combn(posknots, 
            nknots), max(covdata$index))
        notAble <- not_able_knots(allposknots, delta)
        allposknots = allposknots[, -notAble]
        if (ncol(allposknots) == 0) {
            print("allposknots with 0 columns: decrease nknots or delta, or aumentgs sample size\n")
            return(list(knots = NA, min.eqm = NA, model = NULL))
        }
        EQM = NA
        ret = mclapply(1:ncol(allposknots), EQM_cal, bspline, 
            covdata, allposknots, nknots, degree, b.knots, mc.cores = detectCores())
        EQM = simplify2array(ret)
        chosenknots = allposknots[2:(nknots + 1), which.min(EQM)]
        if (isTRUE(bspline)) 
            bestfit = lm(cases ~ bs(index, knots = c(chosenknots), 
                degree = degree), data = covdata)
        else bestfit = lm(cases ~ ns(index, knots = c(chosenknots), 
            Boundary.knots = b.knots), data = covdata)
    }
    else {
        chosenknots = NA
        if (isTRUE(bspline)) 
            bestfit = lm(cases ~ bs(index, knots = NULL, degree = degree), 
                data = covdata)
        else bestfit = lm(cases ~ ns(index, knots = NULL, Boundary.knots = b.knots), 
            data = covdata)
        EQM = mean(bestfit$residuals^2)
    }
    return(list(knots = chosenknots, min.eqm = min(EQM), model = bestfit))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
