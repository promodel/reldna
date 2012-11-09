\name{eldna}
\alias{eldna}
\title{Function to calculate electrostatic profiles of DNA for a list of organisms.}
\description{The main function to calculate electrostatic profile of sequances. It uses function lseqcurlib1D() 
and could be applied to more than one sequance at once. There are also possibilities to save result 
in one or a few Rdata-files.}
\usage{eldna(loadnames, savenames, both = TRUE, ring = FALSE, bound = c(50, 
    nchar(s) - 50), ref = bound[1] + 1, ...)}
\arguments{
  \item{loadnames}{names of files containing sequences of organisms in fasta-format}
  \item{savenames}{names of files to save results}
  \item{both}{logical, if TRUE: electrostatic profiles of both strands is calculated}
  \item{ring}{logical, if TRUE: sequance is (!!!!!!!!)}
  \item{bound}{define a fragment of interest}
  \item{ref}{integer, reference position. Center of base pair number ref is located at the origin of the coordinate system.}
  \item{\dots}{arguments passed to lseqcurlib1D()}
}

\value{Structure of a class "list" with data of electrostatic potential distribution or none.}
\references{\url{http://promodel.icb.psn.ru/publab/?q=ru}

Kamzolova S.G., Sivozhelezov V.S., Sorokin A.A., Dzhelyadin T.R., Ivanova N.N., 
Polozov R.V. // J. Biomol. Struct. Dyn. 2000. V. 18(3). P. 325-334.}
\author{Anatoly Sorokin, Evgenia Temlyakova, Timur Dzhelyadin}





