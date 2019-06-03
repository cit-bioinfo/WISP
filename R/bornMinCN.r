bornMinCN <- function(data, cl, markers_cutoff_pval_anovatest = 0.05, markers_pval_anovatest_fdr = TRUE){
    cl = as.factor(as.character(cl))
    names(cl) = colnames(data)
    resAnova <- resanova(data, cl)
    if (markers_pval_anovatest_fdr) colpval = "anova.adjpvalues" else colpval = "anova.pvalues"
    signifAnova = rownames(resAnova)[resAnova[, colpval] < markers_cutoff_pval_anovatest]
    allp1 <- unlist(lapply(levels(cl), function(i) {
        thiscl = (cl == i) * 1
        names(thiscl) = names(cl)
        a <- sapply(signifAnova, function(g){stats::t.test(unlist(data[g,]) ~ thiscl, alternative="greater")$p.value})
        names(a) = signifAnova
        h1Proportion(a)
    }))
    names(allp1) <- levels(cl)
    minp1 <- allp1[which.min(allp1)]
    if(minp1 < markers_cutoff_pval_anovatest){ bornMIN <- 5 ; warning("At least one population has very few markers") } else bornMIN <- 20
    bornMIN
}
