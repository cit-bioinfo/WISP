\name{SMAPcount}
\alias{SMAPcount}
\title{SMAPcount}
\usage{SMAPcount = function(SMAPBAM, GTF, bamfilepattern=".bam",tumorOrgPattern="hs",hostOrgPattern="mmu", isPairedEnd=TRUE, strandSpecific=TRUE, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType=c("gene_id", "exon_id")[1], useMetaFeatures=TRUE, allowMultiOverlap=TRUE, largestOverlap=TRUE, countMultiMappingReads=FALSE, primaryOnly=FALSE, countChimericFragments=TRUE, ignoreDup=FALSE, requireBothEndsMapped=FALSE,checkFragLength=FALSE)}
\arguments{
\item{SMAPBAM}{defines the path to the STAR mapping output directory containing all the BAMs.}
\item{GTF}{the path and name of the reference gtf file}
\item{bamfilepattern}{the pattern of the BAM files. ".bam" by default}
\item{tumorOrgPattern}{the pattern of the tumor identifier. "hs" by default}
\item{hostOrgPattern}{the pattern of the host identifier. "mmu" by default}
\item{isPairedEnd}{is paired end sequencing. Default is TRUE.}
\item{strandSpecific}{is strand Specific sequencing. Default is TRUE.}
\item{isGTFAnnotationFile}={is the annotation provided in .gtf format. TRUE by default}
\item{GTF.featureType}={feature type used to select rows in the GTF annotation which will be used for read summarization. "exon" by default}
\item{GTF.attrType}{attribute type in the GTF annotation which will be used to group features (eg. exons) into meta-features (eg. genes): either "gene_id" or  "exon_id". "gene_id" by default.}
\item{useMetaFeatures}{see Feaure Count for details. Default is TRUE}
\item{allowMultiOverlap}{see Feaure Count for details. Default is TRUE}
\item{largestOverlap}{see Feaure Count for details. Default is TRUE}
\item{countMultiMappingReads}{see Feaure Count for details. Default is FALSE}
\item{primaryOnly}{see Feaure Count for details. Default is FALSE}
\item{countChimericFragments}{see Feaure Count for details. Default is TRUE}
\item{ignoreDup}{see Feaure Count for details. Default is FALSE}
\item{requireBothEndsMapped}{require both ends to be mapped. default is FALSE}
\item{checkFragLength}{see Feaure Count for details. Default is FALSE}


}
\value{
	
	A list with 2 components: FCcountsTumor and FCcountsHost containing the count data matrix for tumor and host respectivly.}
\description{
This function gives two separate count data matrix for the tumor and host respectively from bam files using featurecount function at the gene or exon level.
}
\examples{
d=system.file("extdata","Example_STAR_OUTPUT", package = "SMAP")
counts = SMAPcount(SMAPBAM=d, GTF =paste(d, "combined/combined.gtf", sep="/"))
countsTumor= counts$FCcountsTumor
countsHost= counts$FCcountsHost

}
