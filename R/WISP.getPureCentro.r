WISP.getPureCentro = function(data,cl,pureSamples_filtering = TRUE, nb_markers_selection = c("custom", "optim_kappa")[1],   nb_markers_max_perClass = 50, markers_cutoff_pval_anovatest = 0.05, markers_pval_anovatest_fdr = TRUE, markers_cutoff_auc = 0.8, pureSamples_deltaTopWeights = 0.5, plot_heatmap = TRUE, add_markers = NULL, col_purePop = NULL){
    
    datatot = data
    cltot = cl
    totsample = colnames(data)
    ind2keep = 0
    if(pureSamples_filtering == TRUE){
        
        resW = data.frame(1:ncol(data))
        while(length(ind2keep) != nrow(resW)){
            print("Retrieving markers")
            resTest = testmarkers(data,cl,markers_cutoff_pval_anovatest, markers_pval_anovatest_fdr)
            print("Centroid calculation")
            genescentro = centroestim(data, cl, resTest = resTest, nb_markers_selection=nb_markers_selection, nb_markers_max_perClass = nb_markers_max_perClass, markers_cutoff_auc = markers_cutoff_auc, add_markers = add_markers,markers_cutoff_pval_anovatest = 0.05, markers_pval_anovatest_fdr = TRUE)$centroid
            print("Weight estimation")
            resW = WISP.getWeight(as.matrix(data[rownames(genescentro),]), as.matrix(genescentro))
            names(cl) = colnames(data)
            resW$labelhisto = as.factor(cl[rownames(resW)])
            resW=resW[resW$WARNING == "OK",]
            ind2keep = rownames(resW)[which(resW[,"topWeightedClass"] == resW[,"labelhisto"] & abs(resW$deltaTopWeights)> pureSamples_deltaTopWeights)]
            if(length(setdiff(colnames(data), ind2keep)) == 0){
                print(paste("No sample was removed" ,setdiff(colnames(data), ind2keep)))
            } else if(length(setdiff(colnames(data), ind2keep)) == 1){
                print(paste(length(setdiff(colnames(data), ind2keep)), "sample was removed:" ,setdiff(colnames(data), ind2keep)))
            } else {
                print(paste(paste(length(setdiff(colnames(data), ind2keep)), "samples were removed:", collapse=" "), paste(setdiff(colnames(data), ind2keep), collapse=" "), collapse=" "))
                
            }
            data = data[,ind2keep]
            cl = cl[ind2keep]
        }
        print(paste("Final Centroid calculation (on ", paste(ncol(data), sep=" "), " samples)", sep=""))
        resTest = testmarkers(data,cl)
        genescentroRes = centroestim(data, cl, resTest = resTest, nb_markers_selection=nb_markers_selection, nb_markers_max_perClass = nb_markers_max_perClass, markers_cutoff_auc = markers_cutoff_auc, add_markers = add_markers)
        genescentro = genescentroRes$centroid
    } else {
        resTest = testmarkers(datatot,cltot,markers_cutoff_pval_anovatest, markers_pval_anovatest_fdr)
        genescentroRes = centroestim(data=datatot, cl=cltot, resTest = resTest, nb_markers_selection=nb_markers_selection, nb_markers_max_perClass = nb_markers_max_perClass, markers_cutoff_auc = markers_cutoff_auc, add_markers = add_markers,markers_cutoff_pval_anovatest = 0.05, markers_pval_anovatest_fdr = TRUE)
        genescentro = genescentroRes$centroid
        ind2keep = colnames(datatot)
    }
    indremoved = setdiff(colnames(datatot),ind2keep)
    
    
    if(plot_heatmap == TRUE){
        data.before = datatot[rownames(genescentro), ]
        data.before = data.before-rowMeans(data.before)
        
        data.after = data[rownames(genescentro), ]
        data.after = data.after-rowMeans(data.after)
        cl = cl[colnames(data.after)]
        
        nbclasses = nlevels(as.factor(cltot))
        
        if(is.null(col_purePop)){
            BaseColors <- setNames(rainbow(nbclasses), levels(as.factor(cltot)))
            
        } else {
            BaseColors = col_purePop
        }
        
        SelectionColors.after = BaseColors[levels(as.factor(as.character(cl)))]
        
        annot1 = data.frame(Classes = sort(cl))
        ha1 = HeatmapAnnotation(df = annot1,
        col = list("Classes" = SelectionColors.after))
        ht1 = ComplexHeatmap::Heatmap(data.after[,names(sort(cl))], circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), cluster_columns = FALSE,show_column_names = T, cluster_rows = FALSE,row_names_side="left",column_names_side="top", name = "Expression",top_annotation = ha1, row_names_gp = grid::gpar(fontsize = 6),column_names_gp = grid::gpar(fontsize = 7), column_title = "Samples kept")
        
        genesclassList = genescentroRes$marksList
        genesclass = unlist(genesclassList)
        names(genesclass) = rep(names(genesclassList), lapply(genesclassList, length))
        if(sum(duplicated(genesclass))>0){
            genesdup = genesclass[ which(duplicated(genesclass))]
            genesclassuniq = genesclass[-which(duplicated(genesclass))]
            names(genesclassuniq)[which(genesclassuniq %in% genesdup)]="notspecific"
            genesclass=genesclassuniq
        }
        annotgenes = data.frame("Classes_marker" = names(genesclass)[match(rownames(data.before),genesclass)])
        hr = rowAnnotation(df = annotgenes, col = list("Classes_marker" = c(SelectionColors.after,"notspecific"="grey")),width = unit(1, "cm"))
        
        if(length(indremoved)!=0){
            if(length(indremoved)==1){
                data.removed = data.frame(data.before[,indremoved])
                rownames(data.removed) = rownames(data.before)
                colnames(data.removed) =indremoved
                cl.removed = cltot[indremoved]
                SelectionColors.removed = BaseColors[levels(as.factor(as.character(cl.removed)))]
                annot2 = data.frame(Classes = cl.removed)
                ha2 = HeatmapAnnotation(df = annot2,
                col = list("Classes" =  SelectionColors.removed))
                ht2 = ComplexHeatmap::Heatmap(data.removed, circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), cluster_columns = FALSE,show_row_names = FALSE,show_column_names = T, cluster_rows = FALSE,column_names_side="top",top_annotation = ha2, column_names_gp = grid::gpar(fontsize = 7), show_heatmap_legend = FALSE,column_title = "Samples removed")
            } else {
                data.removed = data.before[,indremoved]
                cl.removed = cltot[indremoved]
                SelectionColors.removed = BaseColors[levels(as.factor(as.character(cl.removed)))]
                annot2 = data.frame(Classes = sort(cl.removed))
                ha2 = HeatmapAnnotation(df = annot2,
                col = list("Classes" =  SelectionColors.removed))
                ht2 = ComplexHeatmap::Heatmap(data.removed[,names(sort(cl.removed))], circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), cluster_columns = FALSE,show_row_names = FALSE,show_column_names = T, cluster_rows = FALSE,column_names_side="top",top_annotation = ha2, column_names_gp = grid::gpar(fontsize = 7), show_heatmap_legend = FALSE,column_title = "Samples removed")
            }
            
            ComplexHeatmap::draw(hr+ht1+ht2,gap = unit(1.5, "cm"))
            names(indremoved) = cltot[indremoved]
        } else {
            ComplexHeatmap::draw(hr+ht1)
            
        }
    }
    
    names(ind2keep) = cltot[ind2keep]
    if(nb_markers_selection == "optim_kappa"){
        message("Creation of the file \'table_with_allConditionNumber.txt\'")
    }
    return(list("genescentro"  = genescentro, "indkept"=ind2keep, "indremoved"=indremoved, "genesclasses" = genescentroRes$marksList))
    
}
