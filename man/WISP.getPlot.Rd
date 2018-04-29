\name{WISP.getPlot}
\alias{WISP.getPlot}
\title{WISP.getPlot}
\description{Graphical representations of the sample contingent estimations.}
\usage{WISP.getPlot(w, weight_filtered = FALSE, col_purePop = NULL,
annotsup = NULL, col_annotsup = NULL, plot_type = c("heatmap", "barplot")[1], heatmap_orderedby=NULL, barplot_split = c("by.topweight", "by.annotsup")[1]) }
\arguments{
\item{w}{data.frame of weights: object returned by WISP.getWeight}
\item{weight_filtered}{whether use filtered weights. TRUE by default.}
\item{plot_type}{Type of weights representation : heatmap or barplot.}
\item{col_purePop}{(optional) vector of colors for the pure populations. vector names must refer to pure population names that must be ordered as in w.}
\item{annotsup}{(optional) vector of a supplementary annotation to be displayed under each individual sample. supplementary annotation has to be a qualitative variable. vector names must refer to sample names.}
\item{col_annotsup}{(optional) vector of colors for the supplementary annot levels. vector names must refer to the annotation levels.}
\item{heatmap_orderedby}{(optional) the name of a pure population if the user wants the heatmap displayed in a specific order. If not specified, samples are ordered according to the first column in w.}
\item{barplot_split}{whether the barplot should be splitted by top weight pure populations or by the supplementary annotation levels.}

}
\value{
\item{}{Graphical representations of the sample contingent estimations.}
}
\author{Yuna Blum and Aurelie Kamoun}

\keyword{methods}