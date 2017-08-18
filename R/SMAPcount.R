SMAPcount = function(SMAPBAM, GTF, bamfilepattern=".bam",tumorOrgPattern="hs",hostOrgPattern="mmu",
isPairedEnd=TRUE,
strandSpecific=TRUE,
isGTFAnnotationFile=TRUE,
GTF.featureType="exon",
GTF.attrType=c("gene_id", "exon_id")[1],
useMetaFeatures=TRUE,
allowMultiOverlap=TRUE,
largestOverlap=TRUE,
countMultiMappingReads=FALSE,
primaryOnly=FALSE,
countChimericFragments=TRUE,
ignoreDup=FALSE,
requireBothEndsMapped=FALSE,
checkFragLength=FALSE,
...){
    mappedfiles=paste(SMAPBAM,"/",list.files(SMAPBAM,pattern=bamfilepattern,recursive=T),sep="")
    FC=Rsubread::featureCounts(mappedfiles,
    annot.ext=GTF,
    isPairedEnd=isPairedEnd,
    strandSpecific=strandSpecific,
    isGTFAnnotationFile=isGTFAnnotationFile,
    GTF.featureType=GTF.featureType,
    GTF.attrType=GTF.attrType,
    useMetaFeatures=useMetaFeatures,
    allowMultiOverlap=allowMultiOverlap,
    largestOverlap=largestOverlap,
    countMultiMappingReads=countMultiMappingReads,
    primaryOnly=primaryOnly,
    countChimericFragments=countChimericFragments,
    ignoreDup=ignoreDup,
    requireBothEndsMapped=requireBothEndsMapped,
    checkFragLength=checkFragLength,
    ...)
    FCcounts = FC$counts
    hostID = FC$annotation$GeneID[grep(hostOrgPattern, FC$annotation$Chr)]
    tumorID = FC$annotation$GeneID[grep(tumorOrgPattern, FC$annotation$Chr)]
    FCcountsTumor = FCcounts[tumorID,]
    FCcountsHost = FCcounts[hostID,]
    
    save(FCcountsTumor, file="FCcountsTumor.RData")
    save(FCcountsHost, file="FCcountsHost.RData")
    
    return(FCcounts = list(FCcountsTumor=FCcountsTumor, FCcountsHost=FCcountsHost))
    
}
