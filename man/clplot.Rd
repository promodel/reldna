\name{clplot}
\alias{clplot}
\title{Function to plot loadings of components}

\usage{clplot(resPLS, comps = 1:5, legend = "bottomright", lwd = 2, 
    abline = TRUE, save = FALSE, filetype, filename = "PLS_loadings", 
    show = TRUE, ...)}
\arguments{
  \item{resPLS}{mvr-object}
  \item{comps}{numeric vector, numbers of components to show on the plot}
  \item{legend}{character string, specifying the position of the legend}
  \item{lwd}{integer, width of lines}
  \item{abline}{logical, line for TSS-position}
  \item{save}{logical, to save plot into file or not?}
  \item{filetype}{if save=TRUE, specifies file type for saving; possible formats are png, pdf, postscript, jpeg, bmp.}
  \item{filename}{character string, if save=TRUE specifies a name of the file}
  \item{show}{logical, to show plot or not?}
  \item{\dots}{arguments passed to plot}
}



\author{Anatoly Sorokin, Evgenia Temlyakova, Timur Dzhelyadin}





