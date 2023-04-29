getARvsDiff <- function( tmpMatI, tmpMat , n_row, n_col){
  tmpAR <- vector()
  tmpLatDiff <- vector()
  tmpMatrix<- matrix( 0 , n_row, n_col )
  
  for( i in 1:dim( tmpMat)[1]){
    if( is.na( tmpMat[ i, 2])){
      if( length( tmpAR ) == 0 ){
        tmpAR[ length( tmpAR ) + 1 ] <- tmpMat[i,1]
      }
      else{
        print( length( tmpLatDiff ))
        tmpMatrix[ length( tmpAR ) , ] <- tmpLatDiff 
        print( max( abs( tmpLatDiff )))
        tmpAR[ length( tmpAR ) + 1 ] <- tmpMat[i,1]
        tmpLatDiff <- vector()
      }
    }
    else{
      tmpLatDiff[ length( tmpLatDiff) + 1 ] <- car2sph( tmpMat[ i , 1] , tmpMat[ i , 2 ] , tmpMat[ i , 3] )[1,2] - car2sph( tmpMatI[ i , 1] , tmpMatI[ i , 2 ] , tmpMatI[ i , 3] )[1,2]
    }
  }
  tmpMatrix[ length( tmpAR ) , ] <- tmpLatDiff 
  tmpMatrix <- cbind( tmpAR , tmpMatrix )
  tmpLatDiff <- vector()
  return( abs(tmpMatrix) )
}

ce_tmpi <- read.csv("./out_MTPI_cel_MTVariable.csv" , header = F)
ce_tmp <- tmpi <- read.csv("./out_MTP_cel_MTVariable.csv" , header = F)

ce_arvsdiff <- getARvsDiff( ce_tmpi , ce_tmp, 25, 194)

pdf("./ce_diffAngle.pdf")
plot( ce_arvsdiff[,1] , apply( ce_arvsdiff[ , -1] , 1 , max) , xlim = c(1.0,3.0) , ylim = c(0,90) , xlab = "Aspect Ratio" , ylab = "Max difference of the angle between initial and steady state[degree]" , type = "l")
dev.off()

su_tmpi <- read.csv("./out_MTPI_su_MTVariable.csv" , header = F)
su_tmp <- tmpi <- read.csv("./out_MTP_su_MTVariable.csv" , header = F)

su_arvsdiff <- getARvsDiff( su_tmpi , su_tmp, 7 , 3108)

pdf("./su_diffAngle.pdf")
plot( su_arvsdiff[,1] , apply( su_arvsdiff[ , -1] , 1 , max) , xlim = c(1.0,3.0) , ylim = c(0,90) , xlab = "Aspect Ratio" , ylab = "Max difference of the angle between initial and steady state[degree]", type = "l")
dev.off()
