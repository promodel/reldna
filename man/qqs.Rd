\name{qqs}
\alias{qqs}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Template profile of electrostatic potential around individual base pare. 
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Three-dimensional array of potential values at specific points on the cylinder around base pair. The raidus of the cylinder is 15 Angstrom (A) and the length is 210 A, so it cover the range from -105 to +104 A on Z axis.
}
\usage{data(qqs)}
\format{
  The format is:
 num [1:360, 1:210, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
1.	Polozov, R. V. et al. Electrostatic potentials of DNA. Comparative analysis of promoter and nonpromoter nucleotide sequences. J Biomol Struct Dyn 16, 1135-1143 (1999).
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(qqs)
## maybe str(qqs) ; plot(qqs) ...
}
\keyword{datasets}
