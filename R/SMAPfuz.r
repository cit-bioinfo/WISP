#use Rscript --vanilla /datacit/Development/RNAseq_YR/SMAP_fusion/fusionAnalysis_forpaper.r $STAR_FUSION_DIR $SAMPLE
#exemple /datacit/Development/RNAseq_YR/SMAP_fusion/fusionAnalysis_forpaper.r /datacit/Development/RNAseq_YR/SMAP_fusion/ Sample_3

#STAR_FUSION_DIR=commandArgs(trailingOnly=TRUE)[1]
#SAMPLE=commandArgs(trailingOnly=TRUE)[2]

.fuzprocess=function(f,sample){
    
    y=utils::read.delim(f,sep="\t",header=F,as.is=T,comment.char="",skip=1)
    colnames(y)=gsub("^#","",unlist(strsplit(scan(f,sep='\n',what="character",nlines=1,quiet = T)[1],"\t")))
    
    suppressWarnings({
        y[,2]=as.integer(	y[,2])
        y[,3]=as.integer(	y[,3])
    })
    
    y=y[which(!is.na(y[,2])&!is.na(y[,3])),]
    y[,1]=gsub("^#","" 	,y[,1])
    
    y=y[,which(colSums(is.na(y))<nrow(y) )]
    genes=do.call(rbind,strsplit(y[,1],"--"))
    y = data.frame("LeftGene"=genes[,1],   "RightGene"=genes[,2] ,y,stringsAsFactors=F)
    
    y$OrgLeftGene = stringr::str_extract(y$LeftBreakpoint,"mmu|hs")
    y$OrgRightGene = stringr::str_extract(y$RightBreakpoint,"mmu|hs")
    
    y$sample = rep(sample,nrow(y))
    y=unique(y)
    y[order(y[,2],decreasing=T),]
}

.fishersMethod = function(x) stats::pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
# Function to filter dupplicated fusion (selmaxreads)
.dupFUZ=function(expdat){
    expdat$orderfuz =apply(expdat,1,function(x) paste(as.character(sort(x[c("LeftGene","RightGene")])),collapse="--"))
    expdat$sum=apply(expdat[,c("JunctionReads", "SpanningFrags")],1,sum)
    sym=expdat$orderfuz
    i=order(expdat$sum,decreasing=T)
    oki=i[which(!duplicated(sym[i]))]
    newexpdat=data.frame(expdat[oki,])
    rownames(newexpdat)=sym[oki]
    newexpdat
}

#  utils
.getlast=function(x){return(x[length(x)])}
.cmess=function(...,sep=" "){message(paste(...,sep=sep))}
.cstop=function(...,sep=" "){stop(paste(...,sep=sep))}


#' Title
#'  SMAPfuz
#'
#'
#'
#' @param STAR_FUSION_DIR defines the path to the STAR fusion output directory.
#' @param specieOfInterestPattern defines for which specie the fusion are analysed.
#' usually human for Paitent Derived Xenograft
#'
#' @description This function aims at defining the best thresholds for cimeric transcripts detection.
#' It is based on an estimation of an H0 distribution of parameters of eahc of the detected fusion transcritps by
#' taking into account the impossibility of observing cross-species fusions and therefore using them
#' to define the parameters of H0.
#' As of now, only the number of spanning fragments and of junction reads are used.
#'
#'
#' @return A data frame containing one fusion per line.
#' In addition to the information given by STAR fusion, 4 colums are added :
#' JunctionReads_Pvalue : the p-value of the junction read
#' SpanningFrags_Pvalue : the p-value of the spanning fragments
#' Combined_pvalue : The combined p-value using fishers method
#' Combined_adj.pvalue : The FDR correction of the combined p-value
#'
#' @examples
#' d=system.file( "extdata","Example_StarFusion_OUTPUT", package = "SMAP")
#'
#' fusions=SMAPfuz(STAR_FUSION_DIR=d,"hs")
#' fusions_filtered=fusions[which(fusions$Combined_adj.pvalue < 0.01),]
#' fusions_filtered=fusions[which(fusions$Combined_adj.pvalue < 0.001),]
#'



SMAPfuz = function(STAR_FUSION_DIR,specieOfInterestPattern="hs"){
    
    
    finalfusionFiles="fusion.fusion_candidates.final.abridged"
    unfiltfusionFiles="fusion.fusion_candidates.preliminary.filt_info.abridged"
    
    sample =  .getlast(strsplit(STAR_FUSION_DIR,"/")[[1]])
    finaldatafile=grep(finalfusionFiles,list.files(STAR_FUSION_DIR),value=T)
    unfiltdatafile=grep(unfiltfusionFiles,list.files(STAR_FUSION_DIR),value=T)
    
    if(length(finaldatafile)==0 | length(unfiltdatafile)==0 ){
        .cstop("No STAR FUSION files found in the given directory: ",STAR_FUSION_DIR)
    }
    
    
    FilteredFuz=.dupFUZ(.fuzprocess(file.path(STAR_FUSION_DIR,finaldatafile),sample))
    UnfiltFuz=.dupFUZ(.fuzprocess(file.path(STAR_FUSION_DIR,unfiltdatafile),sample))
    
    
    # HO distribution
    valuecols=c("JunctionReads", "SpanningFrags")
    leftgeneOrgCol="OrgLeftGene"
    rightgeneOrgCol="OrgRightGene"
    
    
    
    humanfuze=FilteredFuz[which(FilteredFuz[,rightgeneOrgCol]==specieOfInterestPattern & FilteredFuz[,leftgeneOrgCol]==specieOfInterestPattern),]
    
    reads_junc = UnfiltFuz[which(UnfiltFuz[,rightgeneOrgCol] != UnfiltFuz[,leftgeneOrgCol]),valuecols[1]]
    frag_span = UnfiltFuz[which(UnfiltFuz[,rightgeneOrgCol] != UnfiltFuz[,leftgeneOrgCol]),valuecols[2]]
    nb0junc = sum(reads_junc>0)
    nb0span = sum(frag_span>0)
    #print(nb0junc)
    #print(nb0span)
    
    if(nb0junc<4 & nb0span<4) .cstop("Too few reads support the interspecies-Fusion for H0 calculation")
    
    if(nb0junc> 2 & nb0span<4){
        print("H0 on junction reads only")
        nulld1=MASS::fitdistr(UnfiltFuz[which(UnfiltFuz[,rightgeneOrgCol] != UnfiltFuz[,leftgeneOrgCol]),valuecols[1]],"negative binomial")$estimate
        humanfuze[[paste(valuecols[1],"_","Pvalue",sep="")]] = 1-stats::pnbinom(humanfuze[[valuecols[1]]],size=nulld1["size"],mu=nulld1["mu"])
        humanfuze[[paste(valuecols[1],"_","adj.pvalue",sep="")]]=stats::p.adjust( humanfuze[[paste(valuecols[1],"_","Pvalue",sep="")]],method="fdr")
        return(humanfuze)
    }
    if(nb0junc>2 & nb0span>2){
        print("H0 on junction reads and Spanning Fragments")
        nulld1=MASS::fitdistr (UnfiltFuz[which(UnfiltFuz[,rightgeneOrgCol] != UnfiltFuz[,leftgeneOrgCol]),valuecols[1]],"negative binomial")$estimate
        nulld2=MASS::fitdistr (UnfiltFuz[which(UnfiltFuz[,rightgeneOrgCol] != UnfiltFuz[,leftgeneOrgCol]),valuecols[2]],"negative binomial")$estimate
        logx=unique(floor(2^(seq(0,10,by=0.1))))
        x=log2(logx+1)
        xd1 <- stats::dnbinom(logx, mu =nulld1["mu"], size = nulld1["size"])
        xd2 <- stats::dnbinom(logx, mu = nulld2["mu"], size = nulld2["size"])
        humanfuze[[paste(valuecols[1],"_","Pvalue",sep="")]] = 1-stats::pnbinom(humanfuze[[valuecols[1]]],size=nulld1["size"],mu=nulld1["mu"])
        humanfuze[[paste(valuecols[2],"_","Pvalue",sep="")]]=1-stats::pnbinom(humanfuze[[valuecols[2]]],size=nulld2["size"],mu=nulld2["mu"])
        humanfuze$Combined_pvalue=apply(humanfuze[,c(paste(valuecols[1],"_","Pvalue",sep=""),paste(valuecols[2],"_","Pvalue",sep=""))],1,.fishersMethod)
        humanfuze$Combined_adj.pvalue=stats::p.adjust(humanfuze$Combined_pvalue,method="fdr")
        X=humanfuze
        X$decisionPvals=cut(X$Combined_adj.pvalue,breaks=c(0,0.001,0.01,0.05,0.1,1),include.lowest=T)
        FPvals=cbind(UnfiltFuz[which(UnfiltFuz[,rightgeneOrgCol] != UnfiltFuz[,leftgeneOrgCol]),valuecols[1]],
        UnfiltFuz[which(UnfiltFuz[,rightgeneOrgCol] != UnfiltFuz[,leftgeneOrgCol]),valuecols[2]])
        colnames(FPvals)=valuecols
        xmax=ceiling(log2(1+max(UnfiltFuz[,valuecols[1]])))
        logx=unique(floor(2^(seq(0,xmax,by=0.1))))
        xd1 <- stats::dnbinom(logx, mu =nulld1["mu"], size = nulld1["size"])
        xx=log2(logx+1)
        ymax=ceiling(log2(1+max(UnfiltFuz[,valuecols[2]])))
        logy=unique(floor(2^(seq(0,ymax,by=0.1))))
        xd2 <- stats::dnbinom(logy, mu = nulld2["mu"], size = nulld2["size"])
        xy=log2(logy+1)
        #pdf(paste(STAR_FUSION_DIR,"/", SAMPLE,"/SMAPFusions_fusionsAnalysis_filternDecisionThreshPval.pdf", sep=""))
        
        par(fig=c(0,0.8,0,0.75),mar=c(5.1,4.1,1,1))
        smoothScatter(log2(FPvals+1),transformation = function(x) x^.2,xlim=c(0,ymax),ylim=c(0,xmax) )
        points(log2(X[,valuecols]+1),col=rev(c("black","purple","red","orange","green"))[as.integer(X$decisionPvals)],pch=16)
        legend("topright",paste(c(">0.1","<0.1","<0.05","<0.01","<0.001"),"(fdr)"),col=c("black","purple","red","orange","green"),pch=16, cex=0.9)
        par(fig=c(0.8,1,0,0.75),new=T,mar=c(5.1,1,4.1,1))
        barplot(xd1,xx,type="l",axes=T,bty="l",horiz=T,xlab="Prob. distribution")
        par(fig=c(0,0.8,0.75,1),new=T,mar=c(1,4.1,3,1))
        barplot(xd2,xy,type="l",axes=T,bty="l",main="Parametric p-values (negative binomial)",ylab="Prob. distribution")
        #dev.off()
        
    }
    
    
    
    # isfilt=unlist(lapply(1:nrow(humanfuze),function(i){
    #  thisfuz=humanfuze[i,]
    #  js=which(filtfuz$X.fusion_name == thisfuz$X.fusion_name &
    #             filtfuz$sample == thisfuz$sample &filtfuz$LeftBreakpoint == thisfuz$LeftBreakpoint
    #           &filtfuz$RightBreakpoint == thisfuz$RightBreakpoint&filtfuz$Splice_type == thisfuz$Splice_type )
    #  length(js)
    #}))
    # humanfuze$PassStarFusion =isfilt>0
    return(humanfuze)
    
}


