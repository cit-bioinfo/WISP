\name{WISP-package}
\alias{WISP-package}
\alias{WISP}
\docType{package}
\title{
Weighted In Silico Pathology
}
\description{
an R package for retrieving pure population samples and estimating in a sample the proportion of the different contingents based on the pure profiles.}

\details{
\tabular{ll}{
Package: \tab WISP\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2018-04-07\cr
License: \tab GPL>2\cr
}

}

\author{
Yuna Blum, Aurelien de Reynies\cr
Maintainer: Yuna Blum <yuna.blum@gmail.com>
}
% \references{
% ~~ Literature or other references for background information ~~
% }
% ~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
% ~~ the R documentation directory ~~
% \keyword{ package }
% \seealso{
% ~~ Optional links to other man pages, e.g. ~~
% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
% }
\examples{
data(dataWISP)

resPure = WISP.getPureCentro(data= dataWISP$datapure,cl= dataWISP$clpure, nb_markers_selection = "custom",   nb_markers_max_perClass = 50, pureSamples_deltaTopWeights = 0.6, plot_heatmap = TRUE, col_purePop = c("Epure"="blue","normal"="grey","Spure"="green"))

resW = WISP.getWeight(dataWISP$data, resPure$genescentro)

WISP.getPlot(resW, annotsup = dataWISP$histo, plot_type = "heatmap", col_purePop = c("Epure"="blue","normal"="grey","Spure"="green"))
WISP.getPlot(resW, annotsup = dataWISP$histo, plot_type = "barplot", col_purePop = c("Epure"="blue","normal"="grey","Spure"="green"))

}


\keyword{Deconvolution, Gene expression, contingent.}
