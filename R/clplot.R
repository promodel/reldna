clplot <- function #Function to plot loadings of components
### 
###
(resPLS, ##<< mvr-object
comps = 1:5, ##<< numeric vector, numbers of components to show on the plot
legend = "bottomright", ##<< character string, specifying the position of the legend
lwd = 2, ##<< integer, width of lines
abline = TRUE, ##<< logical, line for TSS-position
save = FALSE, ##<< logical, to save plot into file or not?
filetype, ##<< if save=TRUE, specifies file type for saving; possible formats are png, pdf, postscript, jpeg, bmp.
filename='PLS_loadings', ##<< character string, if save=TRUE specifies a name of the file
show = TRUE, ##<< logical, to show plot or not?
... ##<< arguments passed to plot
) {
    
    if (comps[1] == 1) {
        main <- paste("Loadings of the first ", length(comps), " components", collapse = "")
    }
    if (comps[1] != 1) {
        main <- paste("Loadings of ", comps[1], "-", length(comps), " components", 
            collapse = "")
    }
        cl<-function(){
	  plot(resPLS, "loadings", comps = comps, legend = legend, main = main, lwd = lwd, 
            ...)
	  if (abline == TRUE) {                                                          
	      abline(v = 541, col = "darkgrey")                                          
	  }    
	}
    if (show == TRUE) {
        cl()
    }                                                                                  
    if (save == TRUE) {                                                                
        filesave(cl(), filetype, filename)
    }
    
}
