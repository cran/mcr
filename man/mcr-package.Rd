\docType{package}
\name{mcr-package}
\alias{mcr}
\alias{mcr-package}
\title{Method Comparison Regression}
\description{
 Regression methods to quantify the relation between two measurement methods are provided by this package. 
 In particular it addresses regression problems with errors in both variables and without repeated measurements.
 It implements the CLSI recommendations for analytical method comparison and bias estimation using patient samples. 
 
 The main function for performing regression analysis is \code{\link{mcreg}}. Various functions for summarizing
 and plotting regression results are provided (see examples in \code{\link{mcreg}}).

 For user site testing (installation verification) please use the test case suite provided with the package. 
 The test case suite can be run by sourcing the 'runalltests.R' script in the 'unitTests' folder. It requires the XML
 and Runit packages.
}
\details{
  \tabular{ll}{
    Package: \tab mcr\cr 
    Type: \tab Package\cr
    Version: \tab 1.3.0\cr 
    Date: \tab 2022-09-13\cr 
    License:\tab GPL 3\cr 
    LazyLoad: \tab yes\cr 
  }
}
\author{
  Sergej Potapov <sergej.potapov@roche.com>,
  Fabian Model <fabian.model@roche.com>,
  Andre Schuetzenmeister <andre.schuetzenmeister@roche.com>,
  Ekaterina Manuilova <ekaterina.manuilova@roche.com>,
  Florian Dufey <florian.dufey@roche.com>,
  Jakob Raymaekers <jakob.raymaekers@kuleuven.com> 
}
\references{
  CLSI EP09 \url{https://clsi.org/}
}
\keyword{package}

