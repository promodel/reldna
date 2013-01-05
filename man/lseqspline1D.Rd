\name{lseqspline1D}
\alias{lseqspline1D}
\title{lseqspline1D}
\description{spline function to calulate the profile of electrostatics for long sequences
=======================================================}
\usage{lseqspline1D(s, bound, width = 1, ref)}
\arguments{
  \item{s}{DNA sequence}
  \item{bound}{define a fragment of interest}
  \item{width}{smoothing window width}
  \item{ref}{reference position}
}

\value{list with eight components:
\item{mpot}{1D profile of electrostatic potential along Z axis of DNA;} 
\item{risem}{coordinate of the base pair geometrical center on Z axis of DNA;} 
\item{i1}{index of the base pair geometrical center nearest mesh point ;} 
\item{x}{Z coordinates;} 
\item{seq}{DNA sequence used to calculate profile;} 
\item{bound}{boundaries of the part of interest within the sequence;} 
\item{ref}{index of the base pair that suppose to be placed at the origin;} 
\item{zmap}{data frame of geometical properties of base pairs like index, coordinate of the center, part of profile influenced by its charges.} }

\author{Anatoly Sorokin, Evgenia Temlyakova, Timur Dzhelyadin}





