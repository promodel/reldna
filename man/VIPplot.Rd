\name{VIPplot}
\alias{VIPplot}
\title{Function to plot VIP-data}

\usage{VIPplot(VIPdata, ncomp, col = "darkgrey", drawlines = TRUE, cutoff = seq(0.7, 
    1.1, 0.1), cutoffcol = rainbow(length(cutoff)), numlabels = TRUE, 
    main = "VIP plot", ylab = paste("VIP[", ncomp, "]", sep = ""), 
    xlab = "Angstrom", lwd = 2, cex = 1, save = FALSE, filetype, 
    filename = "VIP", show = TRUE, ...)}
\arguments{
  \item{VIPdata}{matrix with VIP-values}
  \item{ncomp}{integer, the number of component to plot VIP-values}
  \item{col}{color of plot body}
  \item{drawlines}{logical, if TRUE: lines of cutoffs are drawn}
  \item{cutoff}{numeric vector, specifies cutoffs to be drawn (optional)}
  \item{cutoffcol}{character vector, specifies colors for cutoffs (optional)}
  \item{numlabels}{logical, if TRUE: numbers of variables higher or equal to cutoffs' values are shown}
  \item{main}{character string, the title of the plot}
  \item{ylab}{character string, the y label}
  \item{xlab}{character string, the x label}
  \item{lwd}{integer, if drawlines=TRUE: specifies width of lines (optional)}
  \item{cex}{integer, if drawlines=TRUE: specifies symbol size (optional)}
  \item{save}{logical, to save plot into file or not?}
  \item{filetype}{if save=TRUE, specifies file type for saving; possible formats are png, pdf, postscript, jpeg, bmp.}
  \item{filename}{character string, if save=TRUE specifies a name of the file}
  \item{show}{logical, to show plot or not?}
  \item{\dots}{arguments passed to plot}
}



\author{Anatoly Sorokin, Evgenia Temlyakova, Timur Dzhelyadin}





