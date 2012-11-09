filesave<-function#Function for saving plots into a few file formats
### It is used in bodies of several plot functions of the package. Possible formats 
### are png, pdf, postscript, bmp, jpeg. There are some special settings for png to make 
### good-for-our-wiki look of saved picture. 
(imagefun, ##<< function, which output is a plot
filetype, ##<< character string, specifies file type for saving; possible formats are png, pdf, postscript, jpeg, bmp.
width=1000, ##<< width of the plot
height=1000 ##<< height of the plot
){
    if (filetype=='png'){
      filename<-'explvar.png'
      png(filename, bg="white", pointsize=16, width=width, height=height)
			imagefun
      dev.off()
    }
    if (filetype=='postscript'){
      filename<-'explvar.ps'
      postscript(filename, width=width, height=height)
			imagefun
      dev.off()
    }
    if (filetype!='png' & filetype!='postscript'){
      filename<-paste('explvar', filetype, sep='.')
      filetype<-get(filetype)
      filetype(filename, width=width, height=height)
			imagefun
      dev.off()
    }
} 
