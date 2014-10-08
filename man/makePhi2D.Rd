\name{makePhi2D}
\alias{makePhi2D}
\title{makePhi2D}
\description{Function to convert results of APBS multivalue function into
2D matrix for further analysis.}
\usage{makePhi2D(file, xC = 0, yC = 0, zC = 0)}
\arguments{
  \item{file}{input file, result of multivalue invocation.}
  \item{xC}{Position of the  cyllinder center (first number in fgcent)}
  \item{yC}{Position of the cyllinder center (second number in fgcent)}
  \item{zC}{Position of the cyllinder center (last number in fgcent)}
}

\value{matrix with angle in rows and Z in columns}

\author{Anatoly Sorokin, Evgenia Temlyakova, Timur Dzhelyadin}





