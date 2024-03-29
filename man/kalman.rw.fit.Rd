% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.r
\name{kalman.rw.fit}
\alias{kalman.rw.fit}
\title{Cumulative negative log-likelihood}
\usage{
kalman.rw.fit(optim.vars, init.mean.a, init.var.a, x, y, Ts)
}
\arguments{
\item{optim.vars}{Vector of initial values for variables to be estimated (initial mean a,
initial variance of a, b, ln.sig.e, ln.sig.w, Ts)}

\item{init.mean.a}{Starting mean for intercept, a}

\item{init.var.a}{Starting variance for intercept, a}

\item{x}{Independent variable in obs. equation}

\item{y}{Dependent variable in obs. equation}

\item{Ts}{Number of years to omit when calculating the concentrated likelihood for the data set.
See Visser and Molenaar (1988).  Default is 1.}
}
\value{
Cumulative negative log-likelihood
}
\description{
Cumulative negative log-likelihood
}
