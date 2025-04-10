% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection.R
\name{getAMEFit}
\alias{getAMEFit}
\title{Get Fit indices for ame models}
\usage{
getAMEFit(mdl)
}
\arguments{
\item{mdl}{output from ame models of class ame}
}
\value{
a vector of porportions
}
\description{
This function takes an ame model output and
and provides the porportion of simulations for which
the observed value is more distant from the model mean than
the simulated vlaues.
}
