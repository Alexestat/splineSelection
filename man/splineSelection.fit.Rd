\name{splineSelection.fit}
\alias{splineSelection.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
splineSelection.fit(covdata, nlimknot, lambda, delta, degree, bspline = TRUE, b.knots = NULL, nome = "nome")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covdata}{
%%     ~~Describe \code{covdata} here~~
}
  \item{nlimknot}{
%%     ~~Describe \code{nlimknot} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
}
  \item{bspline}{
%%     ~~Describe \code{bspline} here~~
}
  \item{b.knots}{
%%     ~~Describe \code{b.knots} here~~
}
  \item{nome}{
%%     ~~Describe \code{nome} here~~
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
function (covdata, nlimknot, lambda, delta, degree, bspline = TRUE, 
    b.knots = NULL, nome = "nome") 
{
    fixknots = matrix(0, nlimknot + 1, nlimknot + 1)
    if (is.null(b.knots)) {
        if (isTRUE(bspline)) {
            first.knot = sort(covdata$index)[round(length(covdata$index) * 
                0.25)]
            last.knot = max(covdata$index)
            b.knots = c(first.knot, last.knot)
            print(b.knots)
        }
        else {
            first.knot = sort(covdata$index)[round(length(covdata$index) * 
                0.1)]
            last.knot = sort(covdata$index)[round(length(covdata$index) * 
                0.95)]
            b.knots = c(first.knot, last.knot)
        }
    }
    EQMpen = vector()
    for (i in 1:(nlimknot + 1)) {
        a = splineSelection.fitfix(covdata, i - 1, delta, degree, 
            b.knots, bspline)
        EQMpen[i] = a$min.eqm + lambda * i
        fixknots[i, 1:length(a$knots)] = a$knots
    }
    min.nknot = which.min(EQMpen)
    if (min.nknot == 1) {
        bestknots = NULL
    }
    else {
        bestknots = fixknots[which.min(EQMpen), 1:(min.nknot - 
            1)]
    }
    if (isTRUE(bspline)) {
        bestfitall = lm(cases ~ bs(index, knots = c(bestknots), 
            degree = degree), data = covdata)
    }
    else {
        bestfitall = lm(cases ~ ns(index, knots = c(bestknots), 
            Boundary.knots = b.knots), data = covdata)
    }
    return(list(bestknots, bestfitall))
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
