\name{dnaGeom}
\alias{dnaGeom}
\title{dnaGeom}
\description{Function to calculate geometry of the DNA double helix.   }
\usage{dnaGeom(s)}
\arguments{
  \item{s}{sequence of the DNA span to be build}
}

\value{object of class 'DNAgeom', list with four slots:
\item{risem}{coordinate of the base pair geometrical center on Z axis of DNA;} 
\item{twist}{that gives the orientation of the X axis of the base pair;} 
\item{nseq}{simple numerical encoding of a DNA sequence by \code{seqinr::s2n} ;} 
\item{l}{length of the sequence.} }
\references{http://chem.rutgers.edu/~xiangjun/3DNA/images/bp_step_hel.gif}
\author{Anatoly Sorokin, Evgenia Temlyakova, Timur Dzhelyadin}





