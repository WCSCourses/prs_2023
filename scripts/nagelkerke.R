NagelkerkeR2 <- function (rr) 
{
    n <- nrow(rr$model)
    R2 <- (1 - exp((rr$dev - rr$null)/n))/(1 - exp(-rr$null/n))
    RVAL <- list(N = n, R2 = R2)
    return(RVAL)
}