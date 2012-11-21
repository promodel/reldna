#  Derived from MixOmics/R/VIP.R
#  Part of the  MixOmics R package, http://cran.r-project.org/web/packages/mixOmics/index.html
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

VIP <-function#Fuction to calculate VIP-values for variables
##references<< MixOmics
(resPLS  ##<< mvr-object
){
  if(!require(pls)){
    stop('Required library "pls" is missing')
  }
    W<-loadings(resPLS)
    H<-resPLS$ncomp
    q<-ncol(resPLS$model[[names(resPLS$model)[1]]])
    p<-ncol(resPLS$model[[names(resPLS$model)[2]]])
    VIP<-matrix(0, nrow = p, ncol = H)
     
    cor2<-cor(resPLS$model[[names(resPLS$model)[1]]], pls::scores(resPLS), use = "pairwise")^2
    cor2<-as.matrix(cor2, nrow = q)
     
    VIP[, 1]<-W[, 1]^2
     
    if (H > 1) {
        for (h in 2:H) {
            if (q == 1) {
                Rd<-cor2[, 1:h] 
                VIP[, h]<-Rd %*% t(W[, 1:h]^2)/sum(Rd)
            }
            else {
                Rd<-apply(cor2[, 1:h], 2, sum)
                VIP[, h]<-Rd %*% t(W[, 1:h]^2)/sum(Rd)
            }
        }
    }
     
    VIP<-sqrt(p * VIP)
    rownames(VIP)<-rownames(W)
    colnames(VIP)<-paste("comp", 1:H)
     
    return(invisible(VIP))
### Matrix with VIP-values
}


