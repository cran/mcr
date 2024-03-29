% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcPBequi.r
\name{mc.PBequi}
\alias{mc.PBequi}
\title{Equivariant Passing-Bablok Regression}
\usage{
mc.PBequi(
  X,
  Y,
  alpha = 0.05,
  slope.measure = c("radian", "tangent"),
  method.reg = c("PBequi", "TS"),
  extended.output = FALSE,
  calcCI = TRUE,
  methodlarge = TRUE
)
}
\arguments{
\item{X}{measurement values of reference method}

\item{Y}{measurement values of test method}

\item{alpha}{numeric value specifying the 100(1-alpha)\% confidence level}

\item{slope.measure}{angular measure of pairwise slopes  (see \code{\link{mcreg}} for details).\cr
\code{"radian"} - for data sets with even sample numbers median slope is calculated as average of two central slope angles.\cr
\code{"tangent"} - for data sets with even sample numbers median slope is calculated as average of two central slopes (tan(angle)).\cr}

\item{method.reg}{"PBequi" equivariant Passing-Bablok regression; "TS" Theil-Sen regression}

\item{extended.output}{boolean. If TRUE, several intermediate results are returned}

\item{calcCI}{boolean. If TRUE, sd of intercept and slope as well as xw are calculated}

\item{methodlarge}{If TRUE (default), quasilinear method is used, if FALSE, quadratic method is used}
}
\value{
a list with elements.
 \item{b0}{intercept.}
 \item{b1}{slope.}
 \item{se.b0}{respective standard error of intercept.}
 \item{se.b1}{respective standard error of slope.}
 \item{xw}{weighted average of reference method values.}
 \item{weight}{dummy values, only returned it extended.output=FALSE.}
 \item{sez}{variance of intercept for fixed slope (extended.output=TRUE, only).}
 \item{vartau}{variance of Kendall's tau (extended.output=TRUE, only).}
 \item{covtx}{covariance of tau and zeta (extended.output=TRUE, only).}
 \item{x0}{"center of gravity" of x (extended.output=TRUE, only).}
 \item{taui}{"Inversion vector; Indicator of influence"}
}
\description{
This is an implementation of the equivariant Passing-Bablok regression.
}
