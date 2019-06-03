\name{WISP.getWeight}
\alias{WISP.getWeight}
\title{WISP.getWeight}
\description{Pure population weight estimations for each sample.}
\usage{WISP.getWeight(data, centro, scaling = c("none", "scale", "center")[1], cutoff_gobalFtest = 0.05, Rsquared_cutoff = 0.5, cutoff_ttest_weights = 0.05, sum_LessThanOne = TRUE) }
\arguments{
\item{data}{data.frame of gene expression with sample in column.}
\item{centro}{data.frame of centroids given by WISP.getPureClass.}
\item{scaling}{type of data and centroids scaling: none, centered or scaled.}
\item{cutoff_gobalFtest}{p-value threshold of the gobal F-test to considere the model as OK or LIMIT for a sample.}
\item{Rsquared_cutoff}{Rsquared cutoff to considere the model as OK or LIMIT for a sample.}
\item{cutoff_ttest_weights}{p-value threshold of the t-test for each weight.}
\item{sum_LessThanOne}{sum of coefficients less that one (TRUE) or equal to one (FALSE).}
}
\value{
\item{output}{dataframe containing the weight of each pure populations, Distance to model, F.test Pvalue, adjusted R squared and for each weight the Confidence Interval and t.test Pvalue. Last column correpond to a WARNING based on the F-test Pvalue and the adjusted R squared cutoff.}
}
\author{Yuna Blum, Aurelien de Reynies and Sylvie Job}
\keyword{methods}