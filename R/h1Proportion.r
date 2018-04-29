h1Proportion <- function (pv, method = "storey", lambda = 0.5)
{
    method = match.arg(method)
    pi1 = 1 - mean(pv > lambda, na.rm = TRUE)/(1 - lambda)
    if (pi1 < 0) {
        warning(paste("estimated pi1 =", round(pi1, digit = 4),
        "set to 0.0"))
        pi1 = 0
    }
    if (pi1 > 1) {
        warning(paste("estimated pi1 =", round(pi1, digit = 4),
        "set to 1.0"))
        pi1 = 1
    }
    return(pi1)
}
