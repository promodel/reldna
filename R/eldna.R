eldna<-function # Function to calculate electrostatic profiles of DNA for a list of organisms.
### The main function to calculate electrostatic profile of sequances. It uses function lseqcurlib1D() 
### and could be applied to more than one sequance at once. There are also possibilities to save result 
### in one or a few Rdata-files.
##references<< \url{http://promodel.icb.psn.ru/publab/?q=ru}
##references<< Kamzolova S.G., Sivozhelezov V.S., Sorokin A.A., Dzhelyadin T.R., Ivanova N.N., 
## Polozov R.V. // J. Biomol. Struct. Dyn. 2000. V. 18(3). P. 325-334.
(loadnames, ##<< names of files containing sequences of organisms in fasta-format
savenames, ##<< names of files to save results
both=TRUE, ##<< logical, if TRUE: electrostatic profiles of both strands is calculated
ring=FALSE, ##<< logical, if TRUE: sequance is (!!!!!!!!)
bound=c(50, nchar(s)-50), ##<< define a fragment of interest 
ref=bound[1]+1, ##<< integer, reference position. Center of base pair number ref is located at the origin of the coordinate system.   
... ##<< arguments passed to lseqcurlib1D()
){

  savetest<-try(length(savenames), TRUE)
  if (class(savetest) == "try-error"){ savenames<-paste(loadnames, 'eldna', sep='.') }

  if (length(savenames)!=1 & length(loadnames)!=length(savenames)){ stop('length(loadnames)!=length(savenames): specify two vectors with equal length') }

  resultlist<-c()
  k<-0
  for (i in 1:length(loadnames)){
    slist<-read.fasta(file=loadnames[i], as.string=TRUE)
    for (j in 1:length(slist)){
      s<-sub(' ', '', slist[[j]][1])
      bound<-c(50, nchar(s)-50)
      if (ring==TRUE){
	      leftfrag<-substr(s, nchar(s)-bound[1]-1, nchar(s))
	      rightfrag<-substr(s, 1, bound[1])
	      s<-paste(leftfrag, s, rightfrag, sep='')
      }
      result<-lseqcurlib1D(s, bound, ref, ...)
      if (both==TRUE){
				result<-list(forward=result)
				s<-c2s(rev(comp(s2c(s))))
				result$reverse<-lseqcurlib1D(s, bound, ref, ...)
      }
      k<-k+1
      resultlist[[k]]<-result
    }
  }

  if (length(savenames)==1){ 
    filename<-paste(savenames, '.Rdata', sep='')
    save(resultlist, file=filename)   
  }
  if (length(savenames)>1){
    for (j in 1:k){
      filename<-paste(savenames[j], '.Rdata', sep='')
      result<-resultlist[[j]]
      save(result, file=filename)
    }
  }
  return(invisible(resultlist))
### Structure of a class "list" with data of electrostatic potential distribution or none.
} 
