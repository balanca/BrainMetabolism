% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/signal_processing.R
\name{roll.funct}
\alias{roll.funct}
\title{apply a function FUN on a rooling windows of a vector}
\usage{
roll.funct(data, width, FUN, size = T, ...)
}
\arguments{
\item{data}{: a numeric vector}

\item{width}{: the size of the rolling window}

\item{FUN}{: the function to apply}

\item{size}{: a logical value indicating if the returned vector have the same length as original data. default is TRUE.}

\item{...}{: additional argument to pass to the FUN}
}
\value{
a vector with FUN result
}
\description{
apply a function FUN on a rooling windows of a vector
}

