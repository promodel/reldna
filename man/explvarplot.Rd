\name{explvarplot}
\alias{explvarplot}
\title{Function to plot explained variation of a PLS-model}

\usage{explvarplot(resPLS, labels = TRUE, comps = 1:5, save = FALSE, 
    filetype, filename = "Explvar", show = TRUE, cum = TRUE, 
    main = "Explained variation of the model", ...)}
\arguments{
  \item{resPLS}{mvr-object}
  \item{labels}{logical, if TRUE: values of explained variation to be drawn}
  \item{comps}{numeric vector, numbers of components to show on the plot}
  \item{save}{logical, to save plot into file or not?}
  \item{filetype}{if save=TRUE, specifies file type for saving; possible formats are png, pdf, postscript, jpeg, bmp.}
  \item{filename}{character string, if save=TRUE specifies a name of the file}
  \item{show}{logical, to show plot or not?}
  \item{cum}{logical, should it be cumulative sum of explained variation or just an each component contribution}
  \item{main}{character string, specifying title of the plot}
  \item{\dots}{
}
}



\author{Anatoly Sorokin, Evgenia Temlyakova, Timur Dzhelyadin}





