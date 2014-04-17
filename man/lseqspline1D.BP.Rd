\name{lseqspline1D.BP}
\alias{lseqspline1D.BP}
\title{lseqspline1D BP}
\description{spline function to calulate the profile of electrostatics for long 
sequences at centres of base pairs
=======================================================}
\usage{lseqspline1D.BP(s, bound, width = 1, ref)}
\arguments{
  \item{s}{DNA sequence}
  \item{bound}{define a fragment of interest}
  \item{width}{smoothing window width}
  \item{ref}{reference position}
}

\value{list with eight components:
\item{mpot}{1D profile of electrostatic potential along Z axis of DNA;} 
\item{risem}{coordinate of the base pair geometrical center on Z axis of DNA;} 
\item{x}{Z coordinates;} 
\item{seq}{DNA sequence used to calculate profile;} 
\item{bound}{boundaries of the part of interest within the sequence;} 
\item{ref}{index of the base pair that suppose to be placed at the origin;} }

\author{Anatoly Sorokin, Evgenia Temlyakova, Timur Dzhelyadin}





