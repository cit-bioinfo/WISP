\name{SMAPfuz}
\alias{SMAPfuz}
\title{SMAPfuz}
\usage{
SMAPfuz(STAR_FUSION_DIR, specieOfInterestPattern = "hs")
}
\arguments{
\item{STAR_FUSION_DIR}{defines the path to the STAR fusion output directory.}

\item{specieOfInterestPattern}{defines for which specie the fusion are analysed.
usually human for Patient Derived Xenograft}
}
\value{
A data frame containing one fusion per line.
In addition to the information given by STAR fusion, 4 colums are added :
JunctionReads_Pvalue : the p-value of the junction read
SpanningFrags_Pvalue : the p-value of the spanning fragments
Combined_pvalue : The combined p-value using fishers method
Combined_adj.pvalue : The FDR correction of the combined p-value
}
\description{
This function aims at defining the best thresholds for chimeric transcripts detection.
It is based on an estimation of an H0 distribution of parameters of eahc of the detected fusion transcritps by
taking into account the impossibility of observing cross-species fusions and therefore using them
to define the parameters of H0.
As of now, only the number of spanning fragments and of junction reads are used.
}
\examples{
d=system.file( "extdata","Example_StarFusion_OUTPUT", package = "SMAP")

fusions=SMAPfuz(STAR_FUSION_DIR=d,"hs")
fusions_filtered=fusions[which(fusions$Combined_adj.pvalue < 0.01),]
fusions_filtered=fusions[which(fusions$Combined_adj.pvalue < 0.001),]

}
