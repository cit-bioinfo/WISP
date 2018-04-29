resanova <- function(data, cl){
    cl = as.factor(as.character(cl))
    m = rowMeans(data, na.rm = T)
    cl.index = split(1:ncol(data), cl)
    nbclasses = length(cl.index)
    MScond = MSerr = 0
    for (i in 1:nbclasses) {
        cli = cl.index[[i]]
        ni = length(cli)
        x = list(means = rowMeans(data[, cli]), vars = matrixStats::rowVars(as.matrix(data[,
        cli])))
        MScond = MScond + ni * ((x$means - m)^2)
        MSerr = MSerr + (ni - 1) * x$vars
    }
    DF1 = (nbclasses - 1)
    MScond = MScond/DF1
    DF2 = ncol(data) - nbclasses
    MSerr = MSerr/DF2
    N = nrow(data)
    f.stat = matrix(MScond/MSerr, nrow = N)
    pval = stats::pf(f.stat, DF1, DF2)
    PVAL = 1 - pval
    names(PVAL) = rownames(data)
    resAnova = data.frame(anova.pvalues = PVAL, anova.adjpvalues = stats::p.adjust(PVAL, method = "fdr"))
    resAnova
}
