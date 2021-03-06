% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/signal_processing.R
\name{polyfit}
\alias{polyfit}
\title{polynomial fit}
\usage{
polyfit(x, y, order = 1)
}
\arguments{
\item{x}{: a numeric vector}

\item{y}{: a numeric vector (x and y must have the same length)}

\item{order}{: a numeric value for the polynomial degree. default is one.}
}
\value{
model formula

Coef: coeficients in ascending order (i.e. C0, C1, C2, \ldots , Cn )

R2: Rsquare
}
\description{
Fit a n\eqn{^{th}}{th} degree polynom
\deqn{(y = C_nX^n + C_{n-1}X^{n-1} + \ldots + C_1X + C_0)}{(y = CnX^n + Cn-1X^n-1 + \ldots + C1X + C0)}
to vectors x and y
}

