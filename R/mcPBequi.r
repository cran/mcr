###############################################################################
##
## mcPBequi.r
##
## Functions for computing the equivariant Passing Bablok   
## estimator. Unlike the classical Passing Bablok estimator, it can also be 
##calculated for arbitrary slopes neq 1. Furthermore, analytical estimates of 
## the sd of the intercept and slope are provided. 
##
## Copyright (C) 2020 Roche Diagnostics GmbH
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

#' Equivariant Passing-Bablok Regression
#'
#' This is an implementation of the equivariant Passing-Bablok regression. 
#'
#' @param X measurement values of reference method
#' @param Y measurement values of test method
#' @param alpha numeric value specifying the 100(1-alpha)\% confidence level
#' @param slope.measure angular measure of pairwise slopes  (see \code{\link{mcreg}} for details).\cr   
#'          \code{"radian"} - for data sets with even sample numbers median slope is calculated as average of two central slope angles.\cr
#'          \code{"tangent"} - for data sets with even sample numbers median slope is calculated as average of two central slopes (tan(angle)).\cr
#' @param method.reg "PBequi" equivariant Passing-Bablok regression; "TS" Theil-Sen regression
#' @param extended.output boolean. If TRUE, several intermediate results are returned
#' @param calcCI boolean. If TRUE, sd of intercept and slope as well as xw are calculated 
#' @param methodlarge If TRUE (default), quasilinear method is used, if FALSE, quadratic method is used
#' @return a list with elements.
#'  \item{b0}{intercept.}
#'  \item{b1}{slope.}
#'  \item{se.b0}{respective standard error of intercept.}
#'  \item{se.b1}{respective standard error of slope.}
#'  \item{xw}{weighted average of reference method values.}
#'  \item{weight}{dummy values, only returned it extended.output=FALSE.}
#'  \item{sez}{variance of intercept for fixed slope (extended.output=TRUE, only).}
#'  \item{vartau}{variance of Kendall's tau (extended.output=TRUE, only).}
#'  \item{covtx}{covariance of tau and zeta (extended.output=TRUE, only).}
#'  \item{x0}{"center of gravity" of x (extended.output=TRUE, only).}
#'  \item{taui}{"Inversion vector; Indicator of influence"} 
#'


mc.PBequi <- function( X, Y, alpha = 0.05, slope.measure = c("radian", "tangent"), method.reg=c("PBequi","TS"),extended.output=FALSE,calcCI=TRUE,methodlarge=TRUE) {
	## Check validity of parameters
	slope.measure <- match.arg(slope.measure)
	method.reg <- match.arg(method.reg)
    stopifnot(length(X) == length(Y))
	stopifnot(is.numeric(X))
	stopifnot(is.numeric(Y))

	
	######################################################
	n <- length(X)
    z0 <- 1.4 # z0 is the z score for McKean-Schrader
    minvt <- 1E-5 # minimal value for the variance of tau
    maxdConf <- 1.0e0-1.0e-8 # maximal value for dConf
    z <- qt(1-alpha/2,n-2)
    ysign <- 1
    if(method.reg=="PBequi"){
        ysign <- sign(R_ktau(X,Y))
    }
    Y <-Y*ysign

    if (!methodlarge) { #use quadratic method
        pi2 <- atan(Inf)
        intercept <- slope <- 0	
         dx <-outer(X,X,calcDiff)
        diag(dx)<-rep(1,n)
        dy <-outer(Y,Y,calcDiff)
        diag(dy)<-rep(0,n)
        if(method.reg=="PBequi"){
            s<-abs(dy/dx)
        }else{
            s<- dy/dx
        }
        if(slope.measure=="radian") s <- atan(s)
        diag(s) <- NA
        sut <- s[upper.tri(s)]
        sut <- sut[!is.na(sut)]
        if(method.reg=="TS") {
            if(slope.measure=="tangent") sut <- sut[!is.infinite(sut)]
            else sut <- sut[(sut > -pi2) & (sut < pi2)]
        }

        slope <- median(sut)
        s[is.na(s)]<-slope
        if(method.reg=="TS") s[is.infinite(s)] <- slope
        if(method.reg=="TS") {
            if(slope.measure=="tangent") s[is.infinite(s)] <- slope
            else s[(s == -pi2) | (s == pi2)] <- slope
        }
#            print(s-slope)
        if(slope.measure=="radian") slope <- tan(slope)
        
        

        # calculate intercept
        Z <- Y-slope*X
        intercept <- median(Z)


        if(calcCI){
        #   calculate CIs of slope
            if(slope.measure=="radian")slope <- atan(slope)
            taui <- apply(sign(s-slope),1,sum)
            tauii <- n*(n-1)
            vartau <- max(minvt,(sum(taui^2)*4-2*tauii)/n/(n-1)/(n-2)/(n-3))
            CI <- quantile(sut,probs=c(max((1-sqrt(vartau)*z)/2,0),min((1+sqrt(vartau)*z)/2,1)),na.rm=T)
            cil <- as.numeric(CI[1])
            ciu <- as.numeric(CI[2])
            if(slope.measure=="radian") {
                slope <-tan(slope)
                cil <- tan(cil)
                ciu <- tan(ciu)
            }


        #   sd of slope
            seSlope <- (ciu-cil)/(2*z) #McKean-Schrader

            #  standard deviation of intercepts for fixed slope
            SZ <- sort(Z)
            sez <- SZ[n+1-round((n+1)/2-z0*sqrt(n/4),digits=0)]-SZ[round((n+1)/2-z0*sqrt(n/4),digits=0)]

            sez <- sez/(2*z0)



            # covariance of zeta and tau
            Z <- Z-intercept
            if(method.reg=="PBequi"){
                Q <- Y+slope*X
            }else{
                Q <- X
            } 
            ii <- order(Z)
            Z <- Z[ii]
            Q <- Q[ii]
            for (i in 1:(length(Z)-1)){
                   if(calcDiff(Z[i],Z[i+1])==0) Z[i] <- Z[i+1] <- (Z[i]+Z[i+1])/2
                }
            n2 <- n*(n-1)
            zneg <- Z[Z<0]
            zpos <- Z[Z>0]
            qneg <- Q[Z<0]
            qpos <- Q[Z>0]
            nneg <- length(zneg) 
            if(nneg >1){
                nneg2 <- nneg*(nneg-1)
                tauneg <- R_ktau(zneg,qneg)
                tauneg <- tauneg*nneg2
            }else{
                tauneg <- 0
            }
            npos <- length(zpos) 
            if(npos >1){
                npos2 <- npos*(npos-1)
                taupos <- R_ktau(zpos,qpos)
                taupos <- taupos*npos2
            }else{
                taupos <- 0
            }
            covtx <- 2*(taupos-tauneg)/n2/sqrt(n*vartau)
            # derivative of intercept wrt slope
            Zl <- Y-(slope-seSlope*z0)*X
            Zu <- Y-(slope+seSlope*z0)*X
            interceptl <- median(Zl)
            interceptu <- median(Zu)
            x0 <- (interceptl-interceptu)/(2*seSlope*z0)
            # weighted average of X
            xw <- x0-sez/seSlope*covtx 
            # sd of intercept
            seIntercept <- sqrt(max(0,sez^2*(1-covtx^2)+seSlope^2*xw^2))
        }else{
            seIntercept <- NA
            seSlope <- NA
            xw <- NA
            sez <- NA
            vartau <- NA
            covtx <- NA
            x0 <- NA
            taui <- NA
        }
    }else{ # use quasilinear method

        # calculate slope by bisection
        if(method.reg=="PBequi"){
				rcpp_PassingBablok <- utils::getFromNamespace("rcpp_PassingBablok", "robslopes")
                nDupPairs <- countDupPairs(X, Y)
                ind1    <- n * (n - 1) / 2 - nDupPairs
                medind0 <-  upOrdStat(n) # upper median (for intercept calculation)
                medind1  <- upOrdStat(ind1) # upper median
            if(ind1%%2==1){
                slope <- rcpp_PassingBablok(X, Y,FALSE, medind0, medind1)[2]
                #if(slope.measure=="radian") slope <- atan(slope)
            }else{
                slope1 <- rcpp_PassingBablok(X, Y,FALSE, medind0, medind1)[2]
                slope2 <- rcpp_PassingBablok(X, Y,FALSE, medind0, medind1-1)[2]
                if(slope.measure=="tangent"){
                    slope <- (slope1 + slope2)/2 
                }else{
                    slope <- tan((atan(slope1)+ atan(slope2))/2)
                }
            }
        }else{

              # Compute correct order statistics
              nbdups <- countDups(X)
              ind1    <- (n * (n - 1)) / 2 - nbdups 
              medind0        <- upOrdStat(n)
              medind1        <- upOrdStat(ind1)
			  rcpp_TheilSen <- utils::getFromNamespace("rcpp_TheilSen", "robslopes")
              if(((n * (n - 1)) / 2 - nbdups)%%2==1){
                 slope  <- rcpp_TheilSen(X, Y, FALSE, medind0, medind1)[2]
              }else{
                 slope1 <- rcpp_TheilSen(X, Y,FALSE, medind0, medind1)[2]
                 slope2 <- rcpp_TheilSen(X, Y,FALSE, medind0, medind1-1)[2]
                 if(slope.measure=="tangent"){
                     slope <- (slope1 + slope2)/2 
                 }else{
                     slope <- tan((atan(slope1)+ atan(slope2))/2)
                 }
             }
        }




        # calculate intercept
        Z <- Y-slope*X
        intercept <- median(Z)
        if(calcCI){
        #   calculate CIs of slope
#            if(slope.measure=="radian")slope <- atan(slope)

            taui <- tauvar2(X,Y,slope,method.reg=method.reg) # inversion vector
            tauii <- n*(n-1) # count of nonzero pairs
            vartau <- max(minvt,(4*sum(taui^2)-2*tauii)/n/(n-1)/(n-2)/(n-3))

            n2 <- n*(n-1)/2
            dConf <- min(maxdConf,z*sqrt(vartau))
            if(method.reg=="PBequi"){
                cil <- robslopes::PassingBablok(X, Y, alpha = (1-dConf)/2, verbose=FALSE)$slope
                ciu <- robslopes::PassingBablok(X, Y, alpha = (1+dConf)/2, verbose=FALSE)$slope
            }else{
                cil <- robslopes::TheilSen(X, Y, alpha = (1-dConf)/2, verbose=FALSE)$slope
                ciu <- robslopes::TheilSen(X, Y, alpha = (1+dConf)/2, verbose=FALSE)$slope
            }

#            if(slope.measure=="radian") {
#                slope <-tan(slope)
#                cil <- tan(cil)
#                ciu <- tan(ciu)
#            }

        #   sd of slope
            seSlope <- abs((ciu-cil)/(2*z)) #McKean-Schrader

            #  standard deviation of intercepts for fixed slope
            SZ <- sort(Z)
            sez <- SZ[n+1-round((n+1)/2-z0*sqrt(n/4),digits=0)]-SZ[round((n+1)/2-z0*sqrt(n/4),digits=0)]
            sez <- sez/(2*z0)

            # covariance of zeta and tau
            Z <- Z-intercept
            if(method.reg=="PBequi"){
                Q <- Y+slope*X
            }else{
                Q <- X
            } 
            ii <- order(Z)
            Z <- Z[ii]
            Q <- Q[ii]
            for (i in 1:(length(Z)-1)){
                if(calcDiff(Z[i],Z[i+1])==0)  Z[i] <- Z[i+1] <- (Z[i]+Z[i+1])/2
                 }
            n2 <- n*(n-1)
            zneg <- Z[Z<0]
            zpos <- Z[Z>0]
            qneg <- Q[Z<0]
            qpos <- Q[Z>0]
            nneg <- length(zneg) 
            if(nneg >1){
                nneg2 <- nneg*(nneg-1)
                tauneg <- R_ktau(zneg,qneg)
                tauneg <- tauneg*nneg2
            }else{
                tauneg <- 0
            }
            npos <- length(zpos) 
            if(npos >1){
                npos2 <- npos*(npos-1)
                taupos <- R_ktau(zpos,qpos)
                taupos <- taupos*npos2
            }else{
                taupos <- 0
            }
            covtx <- 2*(taupos-tauneg)/n2/sqrt(n*vartau)
            # derivative of intercept wrt slope
            Zl <- Y-(slope-seSlope*z0)*X
            Zu <- Y-(slope+seSlope*z0)*X
            interceptl <- median(Zl)
            interceptu <- median(Zu)
            x0 <- (interceptl-interceptu)/(2*seSlope*z0)
            # weighted average of X
            xw <- x0-sez/seSlope*covtx 
            # sd of intercept
            seIntercept <- sqrt(max(0,sez^2*(1-covtx^2)+seSlope^2*xw^2))
        }else{
            seIntercept <- NA
            seSlope <- NA
            xw <- NA
            sez <- NA
            vartau <- NA
            covtx <- NA
            x0 <- NA
            taui <- NA
        }
    }
    if(method.reg=="PBequi"){
        if(ysign==0){
            intercept <- median(Y)
            slope <- slope*ysign
            covtx <- covtx*ysign
        }else{
            intercept <- intercept*ysign
            slope <- slope*ysign
            covtx <- covtx*ysign
        }
    }

	## Return estimates and sd
    if(extended.output==T){
        ret <- list(b0 = intercept, 
                b1 = slope, 
                se.b0 = seIntercept, 
                se.b1 = seSlope, 
                xw = xw,
                sez = sez,
                vartau=vartau,
                covtx=covtx,
                x0 =x0,
                taui = taui
                )
    }else{
        ret <- list(
                b0 = intercept, 
                b1 = slope, 
                se.b0 = seIntercept, 
                se.b1 = seSlope, 
                xw = xw,
                weight = rep(1,n), 
                taui = taui
        )
    }
    return(ret)
}


tauvar2 <- function(X,Y,m,method.reg=c("PBequi","TS")){ 
    if (method.reg=="PBequi"){
        W <- Y+m*X 
    } else{ #TS
        W <- X
    }
    Z <- Y-m*X
    ok <- complete.cases(W,Z)
    if(sum(ok)>0){
        x1 <- W[ok]
        y1 <- Z[ok]
        ii <- order(x1,y1,decreasing=FALSE)
        xps <- x1[ii]
        yps <- y1[ii]
		countIndividualInversions <- utils::getFromNamespace("countIndividualInversions", "robslopes")
        result <- countIndividualInversions(y = yps)
        result[ii, ] <- result
        taui <- rep(NA,length(x1))
        taui[ok] <- result[, 2] - result[, 1]
    }else{
        taui <- NA
    }
}


countDupPairs <- function(X, Y){
  n <- length(X)
  XY <- cbind(X, Y)
  XY <- XY[order(X, Y), ]
  j <- 1
  jj <- 0
  xy <- XY[1, ]
  for (i in 2:n) {
    if (XY[i, 1] == xy[1] & XY[i, 2] == xy[2]) {
      j <- j + 1
    } else {
      if (j > 1) {
        jj <- jj + j * (j - 1)/2
      }
      j <- 1
      xy <- XY[i, ]
    }
    if (i == n) {
      if (j > 1) {
        jj <- jj + j * (j - 1)/2
      }
    }
  }
  return(jj)
}

upOrdStat <- function(n) {
  # this function returns the order statistic corresponding with the upper median
  if (n %% 2 == 1) {
    return((n + 1) / 2)
  }else{
    return(n / 2 + 1)
  }
}

countDups <- function(x) {
    # this function returns the number of slopes with the same x-coordinate
    x.order     <- order(x)
    xs          <- x[x.order]
    dupnb       <- rle(xs)$lengths
    return(sum(choose(dupnb[dupnb > 1], 2)))
}



