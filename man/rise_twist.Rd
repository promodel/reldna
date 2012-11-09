\name{rise_twist}
\alias{rise_twist}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Parameters of the DNA double helix geometry
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Two main parameters that define geometry of strait double helix are rise, which is the  distance between aromatic plains of two consecutive base pairs, and twist, which is the angle between main axis of of two consecutive base pairs.
}
\usage{data(rise_twist)}
\format{
  A data frame with 16 observations on the following 3 variables.
  \describe{
    \item{\code{dilet}}{a character vector}
    \item{\code{rise}}{a numeric vector}
    \item{\code{twist}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
31.	Ponomarenko M.P., Ponomarenko Iu.V., Kel' A.E., Kolchanov N.A., Karas H., Wingender E., Sklenar H. Computer analysis of conformational features of the eukaryotic TATA-box DNA promotors, Mol Biol (Mosk) 31, 733-740 (1997)
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(rise_twist)
## maybe str(rise_twist) ; plot(rise_twist) ...
}
\keyword{datasets}
