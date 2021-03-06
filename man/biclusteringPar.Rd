% Generated by roxygen2 (4.0.2): do not edit by hand
\name{biclusteringPar}
\alias{biclusteringPar}
\title{Given a set of features, will find the most "dominant" bicluster. Works in
parallel using library "multicore." To set the number of cores used, set
options(cores = N).}
\usage{
biclusteringPar(geneDf, nSamples = 100, lam, lam.lwr = 3.5,
  clustOptions = list())
}
\arguments{
\item{geneDf}{data.frame with genes on rows and columns defining conditions}

\item{nSamples}{integer denoting the number of permutations to perform}

\item{lam}{regularization parameter for conditions (maximum number of conditions allows in a bicluster)}
}
\description{
Given a set of features, will find the most "dominant" bicluster. Works in
parallel using library "multicore." To set the number of cores used, set
options(cores = N).
}

