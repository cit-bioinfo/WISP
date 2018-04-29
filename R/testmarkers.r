testmarkers <- function (data, cl, markers_cutoff_pval_anovatest = 0.05, markers_pval_anovatest_fdr = TRUE)
{
    if (any(is.na(rowMeans(data))))
    stop("Error- missing values in the data")
    resAnova = resanova(data,cl)
    if (markers_pval_anovatest_fdr) {
        colpval = "anova.adjpvalues"
    }
    else {
        colpval = "anova.pvalues"
    }
    if (is.na(sum(resAnova[, colpval]))) {
        resAnova[, colpval][is.na(resAnova[, colpval])] = 1
    }
    if (sum(resAnova[, colpval] < markers_cutoff_pval_anovatest) ==
    0) {
        stop("There is no marker")
    }
    else {
        signifAnova = rownames(resAnova)[resAnova[, colpval] <
        markers_cutoff_pval_anovatest]
        cl = as.factor(as.character(cl))
        names(cl) = colnames(data)
        if (is.null(length(signifAnova)))
        warning("There is no marker")
        else {
            allauc = lapply(levels(cl), function(i) {
                thiscl = (cl == i) * 1
                names(thiscl) = names(cl)
                a = unlist(lapply(signifAnova, function(prob) {
                    PresenceAbsence::auc(data.frame(colnames(data),
                    thiscl[colnames(data)], as.numeric(data[prob,
                    ])), st.dev = F)
                }))
                names(a) = signifAnova
                a
            })
            allauc = do.call(cbind, allauc)
            colnames(allauc) = levels(cl)
            FC = sapply(levels(cl), function(i) {
                classsamples = names(cl)[which(cl == i)]
                othersamples = setdiff(colnames(data), classsamples)
                unlist(lapply(signifAnova, function(p) {
                    mean(as.numeric(data[p, classsamples])) - mean(as.numeric(data[p,
                    othersamples]))
                }))
            })
            colnames(FC) = levels(cl)
            rownames(FC) = signifAnova
        }
        resTest = data.frame(resAnova[signifAnova, ], AUC = allauc,
        logFC = FC)
        return(resTest)
    }
}
