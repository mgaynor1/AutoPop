% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mature.surv.calc.R
\name{mature.surv.calc}
\alias{mature.surv.calc}
\title{Mature survival -  Calculates probability of mature survival.}
\usage{
mature.surv.calc(env.ci, as.msurv)
}
\arguments{
\item{env.ci}{Proportion of environmental variance used to define mature survival rate per generation.
Must be an integer greater than or equal to 0 and less than 1.}

\item{as.msurv}{Mean probability of mature survival. Must be a single integer between 0 and 1.}
}
\value{
Probability of mature survival.
}
\description{
Calculates the probability of mature survival based
on a beta binomial distribution.
}
