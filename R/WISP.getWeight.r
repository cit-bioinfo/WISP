WISP.getWeight = function (data, centro, scaling = c("none", "scale", "center")[1], cutoff_gobalFtest = 0.05, Rsquared_cutoff = 0.2, cutoff_ttest_weights = 0.05)
{
    g = intersect(rownames(data), rownames(centro))
    data = data[g, ]
    centro = as.matrix(centro[g, ])
    if(scaling == "scale"){
        data = scale(data, scale = TRUE)
        centro = scale(centro, scale = TRUE)
    }
    if(scaling == "center"){
        data = scale(data, scale = FALSE)
        centro = scale(centro, scale = FALSE)
    }
    nJ = ncol(centro)
    A = centro
    E = rep(1, nJ)
    F = 1
    G = diag(nrow = nJ)
    H = rep(0, nJ)
    output = data.frame(t(apply(data,2, function(y) {
        B = as.matrix(y)
        out = limSolve::lsei(A, B, E, F, G, H)
        coeff = round(out$X, 4)
        names(coeff) = paste("weight", names(coeff), sep=".")
        vBeta = matrix(coeff)
        dSigmaSq = sum((B - A%*%vBeta)^2)/(nrow(A)-ncol(A))
        dTotalSq = sum((B-mean(B))^2)/(nrow(A)-1)
        dModelSq = sum((A%*%vBeta-mean(B))^2)/(ncol(A)-1)
        mVarCovar = try(dSigmaSq*chol2inv(chol(t(A)%*%A)))
        Adjusted.R.squared = round((dTotalSq-dSigmaSq)/dTotalSq,2)
        Ftest = dModelSq/dSigmaSq
        p.Ftest=stats::pf(q=Ftest, df1=ncol(A)-1, df2=nrow(A)-ncol(A), lower.tail=FALSE)
        ng = nrow(centro)
        
        if(!is.character(mVarCovar)){
            vStdErr = sqrt(diag(mVarCovar))
            CI.inf = sapply(1:nJ,function(j){round(coeff[j] - (stats::qt(1-cutoff_ttest_weights/2,ng-nJ)*vStdErr[j]),2)})
            CI.sup = sapply(1:nJ,function(j){round(coeff[j] + (stats::qt(1-cutoff_ttest_weights/2,ng-nJ)*vStdErr[j]),2)})
            tvalue = sapply(1:nJ,function(j){ coeff[j]/vStdErr[j]})
            p.Ttest = sapply(1:nJ,function(j){stats::pt(abs(coeff[j]/vStdErr[j]), df=ng-nJ, lower=FALSE)*2})
        } else{      mVarCovar = NA
            vStdErr = CI.inf = CI.sup = p.Ttest = rep(NA,nJ)
        }
        names(CI.inf) = paste("CI.inf.", colnames(centro),sep="")
        names(CI.sup) = paste("CI.sup.", colnames(centro),sep="")
        names(p.Ttest) = paste("Pvalue.", colnames(centro),sep="")
        coeff.filtered = coeff
        coeff.filtered[p.Ttest>cutoff_ttest_weights]=0
        names(coeff.filtered) = paste(names(coeff),".filtered",sep="")
        dist.Obs.Model = round(sqrt(sum(((centro %*% coeff) -y)^2)),2)
        c(coeff, dist.Obs.Model = dist.Obs.Model,Ftest.pvalue = p.Ftest, Adjusted.R.squared = Adjusted.R.squared,CI.inf,CI.sup,p.Ttest,coeff.filtered)}
    )))
    output$topWeightedClass = colnames(centro)[apply(output[,
    1:nJ], 1, which.max)]
    output$deltaTopWeights = apply(output[, 1:nJ], 1, function(z) abs(diff(sort(z,
    decreasing = TRUE)[1:2])))
    CI = sapply(1:nJ,function(j){
        CIinf = sapply(output[,paste("CI.inf.", colnames(centro)[j], sep="")], function(x) max(0,x))
        CIsup = sapply(output[,paste("CI.sup.", colnames(centro)[j], sep="")], function(x) min(1,x))
        paste("[",CIinf , ", ", CIsup, "]", sep="")})
    colnames(CI) = colnames(centro)
    output$CI = CI
    output$WARNING = as.factor(c("OK", "LIMIT")[(output$Ftest.pvalue>cutoff_gobalFtest | output$Adjusted.R.squared<Rsquared_cutoff)+1])
    output = output[,-c(grep("CI.inf", colnames(output)), grep("CI.sup", colnames(output)))]
    return(output)
}
