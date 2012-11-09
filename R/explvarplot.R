explvarplot<-function# Function to plot explained variation of a PLS-model
###
###
(resPLS,  ##<< mvr-object
labels=TRUE, ##<< logical, if TRUE: values of explained variation to be drawn
comps=1:5, ##<< numeric vector, numbers of components to show on the plot
save=FALSE, ##<< logical, to save plot into file or not?
filetype, ##<< if save=TRUE, specifies file type for saving; possible formats are png, pdf, postscript, jpeg, bmp.
filename='Explvar', ##<< character string, if save=TRUE specifies a name of the file
show=TRUE, ##<< logical, to show plot or not?
cum=TRUE, ##<< logical, should it be cumulative sum of explained variation or just an each component contribution
main='Explained variation of the model', ##<< character string, specifying title of the plot
...
){

  if (show==TRUE){
    bp<-function(){
      if (cum==TRUE){
        T<-cumsum(explvar(resPLS))[comps]
      }
      else {
        T<-explvar(resPLS)[comps]
      }				
      barplot(T, beside=TRUE, axes=TRUE, ylim=c(0, max(T)+5), xlab='Number of the component', ylab='%', names.arg=seq(1, length(comps)), main=main, ...)
      if (labels==TRUE){ 
	b<-barplot(T, beside=T, axes=TRUE, ylim=c(0, max(T)+5), xlab='Number of the component', ylab='%', names.arg=seq(1, length(comps)), main=main, ...)
	text(b, T+2, round(T, digits = 2)) 
      }
    }
    bp()
  }
	if (save==TRUE) {filesave(bp(), filetype, filename)}
}
