bendability<-function(
  ### spline function to calulate the profile of DNA 
  ### propensity to bend for long sequences according to
  ### I. Brukner, R. Sánchez, D. Suck, S. Pongor, The EMBO Journal 14, 1812–1818 (1995).
  ### =======================================================
  s, ##<< DNA sequence
  bound, ##<< define a fragment of interest
  width=1 ##<< smoothing window width
){
  if (!require(zoo)) {
    stop('Required library "zoo" is not installed.')
  }
  oligos<-getOligos(s,3)
  l<-length(oligos)
  if(length(bound)>2){
    warning(paste('Length of "bound" is',length(bound),'when 2 is expected. First two values of bound are used.'))
    bound<-bound[1:2]
  }else if(length(bound)==1){
    bound<-c(bound,l-bound)
  }
  data<-brukner
  oligos.levels<-tolower(data$Step)
  trin<-factor(oligos, levels = oligos.levels)
  tri<-array(0, dim=c(64, l))
  for(m in 1:l){
    tri[trin[m],m]<-1
  }
  bend<-rollsum(t(data$lnA%*%tri),4,align = 'left')
  return(bend[bound[1]:bound[2]])
}

