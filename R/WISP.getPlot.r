WISP.getPlot <- function(w, weight_filtered = FALSE, col_purePop = NULL,
annotsup = NULL, col_annotsup = NULL, plot_type = c("heatmap", "barplot")[1], heatmap_orderedby=NULL, barplot_split = c("by.topweight", "by.annotsup")[1]){
    
    cl.weight = gsub("weight[.]","", grep("filtered", grep("weight", colnames(w), value = TRUE), value = TRUE, invert = TRUE))
    
    if(is.null(col_purePop)) {
        col_purePop <- setNames(rainbow(length(cl.weight)), cl.weight)
    }
    if(!is.null(annotsup)){
        annotsup.levels = unique(annotsup)
        if(is.null(col_annotsup)) {
            col_annotsup <- setNames(rainbow(length(annotsup.levels)), annotsup.levels)
        }
    }
    if(plot_type=="barplot"){
        if(barplot_split=="by.annotsup" & is.null(annotsup)){
            stop("Error- if barplot_split==by.annotsup, annotsup can't be NULL")
        }
        
        if(weight_filtered == T) {column.weight <- grep("filtered", grep("weight", colnames(w), value = TRUE), value = TRUE)
        } else {column.weight <- grep("filtered", grep("weight", colnames(w), value = TRUE), value = TRUE, invert = TRUE)}
        
        if(barplot_split=="by.topweight"){
            col.class =col_purePop
            annot.samples = setNames(w$topWeightedClass, rownames(w))
            dat <- data.frame(w[, c(column.weight)])
            dat$samplename = rownames(w)
            dat$annot = w$topWeightedClass
        }else {
            col.class =col_annotsup
            annot.samples = annotsup[rownames(w)]
            dat <- data.frame(w[, c(column.weight)])
            dat$samplename = rownames(w)
            dat$annot = annotsup[rownames(w)]
        }
        
        colnames(dat) = gsub("[.]filtered","", gsub("weight[.]","",colnames(dat)))
        suppressMessages({datf <- reshape2::melt(dat)})
        g1=ggplot2::ggplot(data = datf, aes(x = as.factor(samplename), y = value)) +  ylab("Estimated contingent proportions") +
        geom_bar(aes(fill = as.factor(variable)), stat = "identity") +
        facet_grid(~annot, scale="free_x", space="free")+ theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank())+ scale_fill_manual(values=col_purePop)
        
        if(is.null(annotsup)){
            print(g1)
        } else {
            dtannot = data.frame(value = rep(1, length(dat$samplename)))
            dtannot$samplename = dat$samplename
            dtannot$annot = dat$annot
            dtannot$annotsup = annotsup[dat$samplename]
            
            g2=ggplot2::ggplot(data = dtannot, aes(x = as.factor(samplename), y = value)) +
            geom_bar(aes(fill = as.factor(annotsup)), stat = "identity") +
            facet_grid(~annot, scale="free_x", space="free") + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(color="white"),
            axis.text.y=element_text(color="white"),
            axis.ticks.y=element_line(color="white"),
            axis.line = element_line(colour = "white"),
            legend.title=element_blank())+ scale_fill_manual(values=col_annotsup)
            
            print(ggpubr::ggarrange(g1, g2,nrow=2, labels = c("", ""), common.legend = FALSE,heights=c(6, 2),legend="bottom"))
        }
        
    }
    if(plot_type=="heatmap"){
        column.weight <- grep("filtered", grep("weight", colnames(w), value = TRUE), value = TRUE, invert = TRUE)
        if(weight_filtered == T) column.weight <- grep("filtered", grep("weight", colnames(w), value = TRUE), value = TRUE)
        
        if(is.null(heatmap_orderedby)){
            heatmap_orderedby = 1
        } else {
            heatmap_orderedby = grep(heatmap_orderedby,column.weight)
        }
        column.weightdf = data.frame("column.weight" = column.weight)
        
        annotb = data.frame("WARNING" =w[order(w[, column.weight[heatmap_orderedby]]),"WARNING"])
        hb = HeatmapAnnotation(df = annotb,
        col = list("WARNING" = setNames(c("#cb181d","#a1d99b"), c("LIMIT", "OK"))))
        
        hr = rowAnnotation(df = column.weightdf, col = list("column.weight" = setNames(col_purePop,column.weight)),width = unit(1, "cm"))
        if(is.null(annotsup)){
            ht1 = ComplexHeatmap::Heatmap(t(w[order(w[, column.weight[heatmap_orderedby]]),column.weight]), circlize::colorRamp2(c(0, 1), c("white", "black")), cluster_columns = FALSE,show_column_names = T, cluster_rows = FALSE,row_names_side="left",column_names_side="top", name = "weight", bottom_annotation = hb )
            ComplexHeatmap::draw(hr+ht1)
        } else {
            annot = data.frame("annotsup" = annotsup[order(annotsup, w[, column.weight[heatmap_orderedby]])])
            ha = HeatmapAnnotation(df = annot,
            col = list("annotsup" = col_annotsup))
            ht1 = ComplexHeatmap::Heatmap(t(w[rownames(annot),column.weight]), circlize::colorRamp2(c(0, 1), c("white", "black")), cluster_columns = FALSE,show_column_names = T, cluster_rows = FALSE,row_names_side="left",column_names_side="top", name = "weight",top_annotation = ha,row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),bottom_annotation = hb)
            ComplexHeatmap::draw(hr+ht1)
        }
    }
    
    
}
