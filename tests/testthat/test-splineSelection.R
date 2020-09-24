library(splineSelection)

g <- function(x)
    x^2

f <- function(x) {
    g(x) + rnorm(1, sd=.5)
}

test_that("splinefit must select 0 as unique knot", {
    set.seed(1234)
    X <- sample(seq(-2, 2, length=1000), 100, replace = F)
    Y <- unlist(lapply(X, f))
    M <- splinefit(X, Y, 1, 0.2, 1, 1)
    expect_equal(M[[1]], 0)
})
