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
  data<-brukner
  oligos.levels<-tolower(data$Step)
  trin<-factor(oligos, levels = oligos.levels)
  tri<-array(0, dim=c(64, l-1))
  for(m in 1:l-1){
    tri[trin[m],m]<-1
  }
  bend<-rollsum(t(data$lnA%*%tri),4)
  return(bend)
}

