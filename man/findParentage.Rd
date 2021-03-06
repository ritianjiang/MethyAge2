% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findParentage.R
\name{findParentage}
\alias{findParentage}
\title{findParentage}
\usage{
findParentage(DF1, DF2, cors = 4, method = "spearman")
}
\arguments{
\item{DF1}{a data.frame contains methylation percentation}

\item{DF2}{a data.frame contains methylation percentation}

\item{cors}{a numeric value of numbers of cores you wanna use. (default:1)}

\item{method}{whether "spearman" or "pearson}
}
\value{
a dataframe contains pvalue,corr,location
}
\description{
This function find the potential heritable sites between
             p and f1 using correlation test.
}
\examples{
data("ExampleParent")
aaa<-findParentage(Father,Son)
aaa
}
