% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CMSPlus.R
\name{CMSPlus}
\alias{CMSPlus}
\title{Title CMSPlus main function}
\usage{
CMSPlus(exp2symbol, plot = TRUE, parallel.sz = 0)
}
\arguments{
\item{exp2symbol}{a dataframe with Gene Expression Profiles data values, samples in columns, genes in rows, rownames corresponding to gene symbols}

\item{plot}{plot a heatmap for gene signature profile (when plot=TRUE)}

\item{parallel.sz}{Number of processors to use when doing the calculations in parallel. If parallel is loaded and this argument is left with its default value (parallel.sz=0) then it will use all available core processors unless we set this argument with a smaller number.}
}
\value{
a phenotype data frame with both CMSPlus labels and CMSClassifier labels using random forest model
}
\description{
Title CMSPlus main function
}
