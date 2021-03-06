% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/savage.R
\name{savage.select}
\alias{savage.select}
\title{Choose the best plan based on Savage using R}
\usage{
savage.select(x)
}
\arguments{
\item{x}{the data matrix which the columns represent the plans and the rows represent the variables considered}
}
\value{
the best solution based on Savage criteria and its compromise standard return value
}
\description{
Choose the best plan based on Savage using R
}
\examples{
\dontrun{
x<-matrix(c(30,20,10,12,15,9,-6,-2,5),3)
savage.select(x)
}
}
