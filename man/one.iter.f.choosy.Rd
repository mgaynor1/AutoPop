% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/one.iter.f.choosy.R
\name{one.iter.f.choosy}
\alias{one.iter.f.choosy}
\title{Simulate a single generation.}
\usage{
one.iter.f.choosy(
  ct.vec,
  aii.vec,
  as.matur,
  as.msurv.e.set,
  d,
  gnum.base,
  s,
  b,
  cc,
  mc,
  density.type = "all",
  mate.lazy = FALSE
)
}
\arguments{
\item{ct.vec}{Population composition at time t.
Indicates the sum of each type at time t, \verb{ct.vec = c(2x_immature, 3x_immature, 4x_immature, 2x_mature, 3x_mature, 4x_mature)}}

\item{aii.vec}{The survival probability of an immature individual for each cytotype. Must be a list of three integers between 0 and 1. For example, \code{aii.vec = c(0.5, 0.3, 0.5)}.}

\item{as.matur}{The probability of maturation from an immature stage to the mature stage for each cytotype. Must be a list of three integers between 0 and 1. For example, \code{as.matur = c(0.5, 0.3, 0.5)}.}

\item{as.msurv.e.set}{The survival probability of a mature individual for each cytotype. Must be a list of three integers between 0 and 1. For example, \code{as.msurv.e.set = c(0.5, 0.3, 0.5)}.
This list is defined within the \code{gen.iter.f.choosy()} function.}

\item{d}{Strength of density dependency on gamete production for each cytotype. Must be a list of three integers between 0 and 1. For example, \code{d = c(0.001, 0.009, 0.001)}.}

\item{gnum.base}{Mean number of gametes per individual per cytotype. Must be a list of three numeric values. For example, \code{gnum.base = c(100, 100, 100)}.}

\item{s}{Selfing rate. Must be a single integer between 0 and 1.}

\item{b}{Proportion of unreduced gamete formed by each diploid and tetraploid individual.
Must be a single integer between 0 and 1.}

\item{cc}{Proportion of 3n gamete formation from each triploid individual.
Must be a single integer between 0 and 1.}

\item{mc}{Strength of mating choice. Must be a single integer between 0 and 1.}

\item{density.type}{Default = "all", this sets the density at time t as all individuals at time t.
Alternatively, "like-cytotype" sets the density at time t for each cytotype based on only the total immature and mature individuals
of that cytotype at time t.}

\item{mate.lazy}{Default = FALSE, this prevents selfing from occurring during outcrossing. However, this increases the computational time by 31x!}
}
\value{
List of 9, with 1:6 representing the number of individuals of each cytotype
and at both stages at time t + 1. Items 7:9 are the number of gametes sampled for
2x, 3x, and 4x individuals at time t.
}
\description{
Defined in detail in \url{Gaynor et al. 2023}. To summarize,
this function is a single time step in a stochastic stage-structured
matrix population dynamics model for diploid, triploid, and autotetraploid
perennial plants with two life-stages (reproductively immature and reproductively mature).
Population composition at time t + 1 is defined by reproduction,
survival, and maturation. Note, lists of three should always have values representing \code{c(diploids, triploid, autotetraploids)}.
}
