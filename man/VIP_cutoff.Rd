\name{VIP_cutoff}
\alias{VIP_cutoff}
\title{Function to determine VIP-values higher or equal to cutoffs}

\usage{VIP_cutoff(VIPdata, ncomp, cutoffs, save = FALSE, filename = paste("VIP_cutoff_", 
    paste(cutoffs, collapse = "_"), ".Rdata", sep = ""))}
\arguments{
  \item{VIPdata}{matrix with VIP-values}
  \item{ncomp}{integer, the number of component to plot VIP-values}
  \item{cutoffs}{numeric vector, specifies cutoffs}
  \item{save}{logical, to save results into Rdata-file or not?}
  \item{filename}{character string, if save=TRUE specifies a name of the file}
}



\author{Anatoly Sorokin, Evgenia Temlyakova, Timur Dzhelyadin}





