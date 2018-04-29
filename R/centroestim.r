centroestim <- function (data, cl, resTest, nb_markers_selection, nb_markers_max_perClass, markers_cutoff_auc, add_markers = NULL,
markers_cutoff_pval_anovatest = 0.05, markers_pval_anovatest_fdr = TRUE)
{
    cl = as.factor(as.character(cl))
    nbclasses = length(levels(cl))
    if (nb_markers_selection == "custom") {
        marksList = lapply(levels(cl), function(aclass) {
            FC = resTest[resTest[, paste("AUC", aclass, sep = ".")] >=
            markers_cutoff_auc, ]
            rownames(FC)[order(FC[, paste("logFC", aclass, sep = ".")],
            decreasing = T)][1:min(nrow(FC), nb_markers_max_perClass)]
        })
        names(marksList) = levels(cl)
        if (!is.null(add_markers)) {
            marksList$add = add_markers
        }
    }
    if (nb_markers_selection == "optim_kappa") {
        bornMIN <- bornMinCN(data, cl, markers_cutoff_pval_anovatest, markers_pval_anovatest_fdr)
        allmarksList = lapply(bornMIN:200, function(x) {
            marksList = lapply(levels(cl), function(aclass) {
                resTest <- resTest[order(resTest[, paste("AUC",
                aclass, sep = ".")], resTest[, paste("logFC",
                aclass, sep = ".")], decreasing = T), ]
                rownames(resTest)[1:x]
            })
            names(marksList) = levels(cl)
            return(marksList)
        })
        allmarks = lapply(allmarksList, function(x) unique(unlist(x)))
        allcondnumber <- lapply(allmarks, function(marks) {
            centroid = data.frame(t(apply(data[marks, ], 1, function(x) tapply(x,
            cl, mean))))
            kappa(centroid)
        })
        sumcn <- data.frame(geneNumber=unlist(lapply(allmarks, length)), conditionNumber=unlist(allcondnumber))
        write.table(sumcn, file="table_with_allConditionNumber.txt", sep="\t", row.names=FALSE)
        bestmarks <- allmarks[[which.min(allcondnumber)]]
        if (is.null(add_markers)) {
            centroid = data.frame(t(apply(data[bestmarks, ],
            1, function(x) tapply(x, cl, mean))))
            marksList = allmarksList[[which.min(allcondnumber)]]
        }
        else {
            add_markers_present = intersect(rownames(data), add_markers)
            if (length(intersect(rownames(data), add_markers)) ==
            0) {
                centroid = data.frame(t(apply(data[bestmarks,
                ], 1, function(x) tapply(x, cl, mean))))
                marksList = allmarksList[[which.min(allcondnumber)]]
            }
            else {
                centroid = data.frame(t(apply(data[unique(c(bestmarks,
                add_markers_present)), ], 1, function(x) tapply(x,
                cl, mean))))
                marksList = allmarksList[[which.min(allcondnumber)]]
                marksList$add = add_markers_present
            }
        }
    }
    centroid = data.frame(t(apply(data[unique(unlist(marksList)),
    ], 1, function(x) tapply(x, cl, mean))))
    return(list(centroid = centroid, marksList = marksList))
}
