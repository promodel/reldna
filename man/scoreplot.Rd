\name{scoreplot}
\alias{scoreplot}
\title{Function to make scoreplots}

\usage{scoreplot(resPLS, comps = c(1, 2), groups, col, ellipse = TRUE, 
    labels = TRUE, pch = 20, save = FALSE, filetype, filename = "PLS_scores", 
    show = TRUE)}
\arguments{
  \item{resPLS}{mvr-object}
  \item{comps}{numeric vector of length equals two, numbers of components for scoreplot}
  \item{groups}{numeric vector or list, containing indices of observation to show in different colors (optional)}
  \item{col}{if groups are specified, it is possible to choose colors for them (optional)}
  \item{ellipse}{logical, if TRUE: ellipse that contains 95\% of the observed data is drawn}
  \item{labels}{logical, if TRUE: observations are labeled}
  \item{pch}{integer, specifies symbol type in plot}
  \item{save}{logical, to save plot or not?}
  \item{filetype}{if save=TRUE, specifies file type for saving; possible formats are png, pdf, postscript, jpeg, bmp.}
  \item{filename}{character string, if save=TRUE specifies a name of the file}
  \item{show}{logical, to show plot or not?}
}



\author{Anatoly Sorokin, Evgenia Temlyakova, Timur Dzhelyadin}





