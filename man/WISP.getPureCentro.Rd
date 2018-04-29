\name{WISP.getPureCentro}
\alias{WISP.getPureCentro}
\title{WISP.getPureCentro}
\description{Get pure population centroid profiles.}
\usage{WISP.getPureCentro(data,cl,pureSamples_filtering = TRUE, nb_markers_selection = c("custom", "optim_kappa")[1],   nb_markers_max_perClass = 50, markers_cutoff_pval_anovatest = 0.05, markers_pval_anovatest_fdr = TRUE, markers_cutoff_auc = 0.8, pureSamples_deltaTopWeights = 0.5, plot_heatmap = TRUE, add_markers = NULL, col_purePop = NULL) }
\arguments{
\item{data}{data.frame of gene expression restricted to samples considered as pure by the user.}
\item{cl}{vector of pure population class for each sample, where samples are in the same order than in data object. vector names must refer to sample names.}
\item{pureSamples_filtering}{should WISP proceeds its pure sample filtering }
\item{nb_markers_selection}{marker selection method : custom means based on nb_markers_max_perClass, optim_kappa will find the optimal number of markers based on the conditional number.}
\item{nb_markers_max_perClass}{maximum number of markers per class (pure polulations).}
\item{markers_cutoff_pval_anovatest}{cutoff used to retrieve significant genes based on the ANOVA test comparing the different pure populations.}
\item{markers_pval_anovatest_fdr}{Should fdr adjusted P-value be used to test gene significance?}
\item{markers_cutoff_auc}{cutoff used for AUC criteria in the marker selection procedure if nb_markers_selection = custom}
\item{pureSamples_deltaTopWeights}{weight difference between the 2 highest weights for a sample.}
\item{plot_heatmap}{whether display heatmap of kept and removed sample expression profiles based on final marker selection.}
\item{add_markers}{(optional) vector of gene names to be added in the centroid signature.}
\item{col_purePop}{(optional) vector of colors for the pure populations. vector names must refer to pure population names.}
}
\value{The ouput is a list with following objects:
\item{genescentro}{centroids of pure populations based on the pure samples.}
\item{indkept}{pure samples kept by the procedure.}
\item{indremoved}{samples not considered as pure that where removed by WISP.}
\item{genesclasses}{list of markers that were used in centroid calculation for each pure populations.}


}
\author{Yuna Blum, Aurelien de Reynies and Sylvie Job}
\keyword{methods}