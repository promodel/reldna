lseqspline1D<-function(
  ### spline function to calulate the profile of electrostatics for long sequences
  ### !!!!! DO NOT USE. FUNCTION IS UNDER DEVELOPMENT !!!!!!!
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
    }else if(length(bound==1)){
      bound<-c(bound,geom$l-bound)
    }
    zout<-(risem[bound[1]]-9):(risem[bound[2]]-9)
  }else{
    stop('no value for "bound" is provided');
  }
  zlib<-dim(qqs)[2]
  lout<-length(zout)
  exZ<- (-lout-zlib/2):(lout+zlib/2-1)
  pad<-rep(0,lout)  
  qqm<-apply(qqs, c(2,3), FUN=mean)
  
  spq<-list(
    A=splinefun(exZ,c(pad,qqm[,1],pad),method='natural'),
    T=splinefun(exZ,c(pad,qqm[,2],pad),method='natural'),
    G=splinefun(exZ,c(pad,qqm[,3],pad),method='natural'),
    C=splinefun(exZ,c(pad,qqm[,4],pad),method='natural'))
  i1<-floor(risem[bound[1]:bound[2]]-risem[ref]+9)
  
#  mz<-matrix(zout,nrow=length(risem),ncol=lout,byrow=TRUE)
#  mr<-matrix(risem,nrow=length(risem),ncol=lout,byrow = FALSE)
#  msp<-mz-mr
  pot<-rep(0,lout)
  for(i in 1:dim(msp)[1]){
    msp<-zout - risem[i]
    pot<-pot+spq[[geom$nseq[i]]](msp)
  }
  elstatlist<-list(mpot=pot, risef=risem, i1=i1,x=zout,seq=s,bound=bound,ref=ref)
  class(elstatlist)<-'elDNA1d'
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
  qqm<-apply(qqs, c(2,3), FUN=mean)

  spq<-list(
    A=splinefun(exZ,c(pad,qqm[,1],pad),method='natural'),
    T=splinefun(exZ,c(pad,qqm[,2],pad),method='natural'),
    G=splinefun(exZ,c(pad,qqm[,3],pad),method='natural'),
    C=splinefun(exZ,c(pad,qqm[,4],pad),method='natural'))
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
}