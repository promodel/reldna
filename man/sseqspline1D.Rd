\name{sseqspline1D}
\alias{sseqspline1D}
\title{sseqspline1D}
\description{spline function to calulate the profile of electrostatics for short sequences}
\usage{sseqspline1D(s, bound = NA, width = 1, ref, zout = -540:179)}
\arguments{
  \item{s}{DNA sequence}
  \item{bound}{define a fragment of interest. NA enforce the use 'zout' value. Not NA value will override 'zout' value.}
  \item{width}{smoothing window width}
  \item{ref}{reference position}
  \item{zout}{exact coordinates in which values of the potential will be calculated. Not NA 'bound' value will make 'zout' overwritten.}
}



\author{Anatoly Sorokin, Evgenia Temlyakova, Timur Dzhelyadin}





