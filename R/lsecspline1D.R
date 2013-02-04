.makeSplines1D<-function(
  ### internal function to prepare spline function of the electrostatic profile for individual base pair type
  ### =======================================================
){
  qqm<-apply(qqs, c(2,3), FUN=mean)
  zlib<-dim(qqs)[2]
  exZ<- (-zlib/2):(zlib/2-1)
  spq<-list(
    A=splinefun(exZ,qqm[,1],method='monoH.FC'),
    T=splinefun(exZ,qqm[,2],method='monoH.FC'),
    G=splinefun(exZ,qqm[,3],method='monoH.FC'),
    C=splinefun(exZ,qqm[,4],method='monoH.FC'))
  return(spq)
}

lseqspline1D<-function(
  ### spline function to calulate the profile of electrostatics for long sequences
  ### =======================================================
  s, ##<< DNA sequence
  bound, ##<< define a fragment of interest
  width=1, ##<< smoothing window width
  ref ##<< reference position 
  #,filename=NA #<< name of the file to save data in. Empty string or NA value means file would not be saved
  #,name ##<< name of the library
){
  geom<-dnaGeom(s)
  risem<-geom$risem-geom$risem[ref]
  if(all(!is.na(bound))&length(bound)>0){
    if(length(bound)>2){
      warning(paste('Length of "bound" is',length(bound),'when 2 is expected. First two values of bound are used.'))
      bound<-bound[1:2]
    }else if(length(bound)==1){
      bound<-c(bound,geom$l-bound)
    }
    zout<-floor(risem[bound[1]]-9):(risem[bound[2]]+9)
  }else{
    stop('no value for "bound" is provided');
  }
  zlib<-dim(qqs)[2]
  lout<-length(zout)
#  exZ<- (-zlib/2):(zlib/2-1)
  pad<-rep(0,lout)  
  
  spq<-.makeSplines1D()

  sgz<-data.frame(pos=1:length(s),risem=risem,
                  minZ=(risem-zlib/2),
                  maxZ=(risem+zlib/2-1))
  sgz$minZ[sgz$minZ<min(zout)]<-min(zout)
  sgz$maxZ[sgz$maxZ>max(zout)]<-max(zout)
  sgz$minI<-floor(sgz$minZ)-min(zout)+1
  sgz$maxI<-ceiling(sgz$maxZ)-min(zout)+1  
  
  ind<-which(sgz$minZ<sgz$maxZ)
  i1<-floor(risem[bound[1]:bound[2]]-risem[ref]+9)
  
  pot<-rep(0,lout)
  for(i in ind){
    indp<-sgz$minI[i]:sgz$maxI[i]
    pot[indp]<-pot[indp]+spq[[geom$nseq[i]]](zout[indp]-sgz$risem[i])
  }
  elstatlist<-list(mpot=pot, risem=risem, i1=i1,x=zout,seq=s,bound=bound,ref=ref,zmap=sgz[ind,])
  class(elstatlist)<-'elDNA1d'
  ##value<< list with eight components:
  ##\item{mpot}{1D profile of electrostatic potential along Z axis of DNA;} 
  ##\item{risem}{coordinate of the base pair geometrical center on Z axis of DNA;} 
  ##\item{i1}{index of the base pair geometrical center nearest mesh point ;} 
  ##\item{x}{Z coordinates;} 
  ##\item{seq}{DNA sequence used to calculate profile;} 
  ##\item{bound}{boundaries of the part of interest within the sequence;} 
  ##\item{ref}{index of the base pair that suppose to be placed at the origin;} 
  ##\item{zmap}{data frame of geometical properties of base pairs like index, coordinate of the center, part of profile influenced by its charges.} 
  return (elstatlist)
}


sseqspline1D<-function(
  ### spline function to calulate the profile of electrostatics for short sequences
  s, ##<< DNA sequence
  ref, ##<< reference position 
  zout=-540:179 ##<< exact coordinates in which values of the potential will be calculated.
  #,filename=NA #<< name of the file to save data in. Empty string or NA value means file would not be saved
  #,name ##<< name of the library
){
  geom<-dnaGeom(s)
  risem<-geom$risem-geom$risem[ref]
  zlib<-dim(qqs)[2]
  lout<-length(zout)
  exZ<- (-lout-zlib/2):(lout+zlib/2-1)
  pad<-rep(0,lout)  
  spq<-.makeSplines1D()
  mz<-matrix(zout,nrow=length(risem),ncol=lout,byrow=TRUE)
  mr<-matrix(risem,nrow=length(risem),ncol=lout,byrow = FALSE)
  msp<-mz-mr
  pot<-rep(0,lout)
  for(i in 1:dim(msp)[1]){
    pot<-pot+spq[[geom$nseq[i]]](msp[i,])
  }
  return(pot)
}

sseqspline1D.BP<-function(  
### spline function to calulate the profile of electrostatics for short sequences in base pairs
  s, ##<< DNA sequence
  ref, ##<< reference position 
  zout=-540:179 ##<< exact coordinates in which values of the potential will be calculated.
  #,filename=NA #<< name of the file to save data in. Empty string or NA value means file would not be saved
  #,name ##<< name of the library
){
  geom<-dnaGeom(s)
  risem<-geom$risem-geom$risem[ref]
  pot<-sseqspline1D(s,ref,zout)
  ind<-which(risem>=min(zout)&risem<=max(zout))
  bp<-ind-ref
  pspline=splinefun(zout,pot,method='natural')
  potbp<-pspline(risem[ind])
  return(cbind(bp,potbp))
  ### potential profile values at the centre of the base pair
}
