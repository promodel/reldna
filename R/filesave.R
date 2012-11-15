filesave<-function#Function for saving plots into a few file formats
### It is used in bodies of several plot functions of the package. Possible formats 
### are png, pdf, postscript, bmp, jpeg. There are some special settings for png to make 
### good-for-our-wiki look of saved picture. 
(imagefun, ##<< function, which output is a plot
filetype, ##<< character string, specifies file type for saving; possible formats are png, pdf, postscript, jpeg, bmp.
filename='explvar',##<< character string, specifies file name
width=1000, ##<< width of the plot
height=1000 ##<< height of the plot
){
    if (filetype=='png'){
      filename<-paste(filename,'.png',sep='')
      png(filename, bg="white", pointsize=16, width=width, height=height)
			imagefun
      dev.off()
    }
    if (filetype=='postscript'){
      filename<-paste(filename,'.ps',sep='')
      postscript(filename, width=width, height=height)
			imagefun
      dev.off()
    }
    if (filetype!='png' & filetype!='postscript'){
      filename<-paste(filename, filetype, sep='.')
      filetype<-get(filetype)
      filetype(filename, width=width, height=height)
			imagefun
      dev.off()
    }
} 
