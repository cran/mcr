\name{mcreg}
\alias{mcreg}
\title{Comparison of Two Measurement Methods Using Regression Analysis}
\usage{
mcreg(x, y = NULL, error.ratio = 1, alpha = 0.05, mref.name = NULL,
  mtest.name = NULL, sample.names = NULL, method.reg = c("PaBa", "LinReg",
  "WLinReg", "Deming", "WDeming", "PaBaLarge"), method.ci = c("bootstrap",
  "jackknife", "analytical", "nestedbootstrap"),
  method.bootstrap.ci = c("quantile", "Student", "BCa", "tBoot"),
  nsamples = 999, nnested = 25, rng.seed = NULL,
  rng.kind = "Mersenne-Twister", iter.max = 30, threshold = 1e-06,
  na.rm = FALSE, NBins = 1e+06)
}
\arguments{
  \item{x}{measurement values of reference method, or two
  column matrix.}

  \item{y}{measurement values of test method.}

  \item{mref.name}{name of reference method (Default
  "Method1").}

  \item{mtest.name}{name of test Method (Default
  "Method2").}

  \item{sample.names}{names of cases (Default "S##").}

  \item{error.ratio}{ratio between squared measurement
  errors of reference and test method, necessary for Deming
  regression (Default 1).}

  \item{alpha}{value specifying the 100(1-alpha)\%
  confidence level for confidence intervals (Default is
  0.05).}

  \item{method.reg}{regression method.  It is possible to
  choose between five regression methods: \code{"LinReg"} -
  ordinary least square regression.\cr \code{"WLinReg"} -
  weighted ordinary least square regression.\cr
  \code{"Deming"} - Deming regression.\cr \code{"WDeming"}
  - weighted Deming regression.\cr \code{"PaBa"} -
  Passing-Bablok regression.\cr \code{"PaBaLarge"} -
  approximative Passing-Bablok regression for large
  datasets, operating on \code{NBins} classes which each
  slope is classified to instead of building the complete
  triangular matrix of all N*N/2 slopes.}

  \item{method.ci}{method of confidence interval
  calculation. The function contains four basic methods for
  calculation of confidence intervals for regression
  coefficients.  \code{"analytical"} - with parametric
  method.\cr \code{"jackknife"} - with leave one out
  resampling.\cr \code{"bootstrap"} - with ordinary
  non-parametric bootstrap resampling.\cr \code{"nested
  bootstrap"} - with ordinary nonparametric bootstrap
  resampling.\cr}

  \item{method.bootstrap.ci}{bootstrap based confidence
  interval estimation method.}

  \item{nsamples}{number of bootstrap samples.}

  \item{nnested}{number of nested bootstrap samples.}

  \item{rng.seed}{integer number that sets the random
  number generator seed for bootstrap sampling. If set to
  NULL currently in the R session used RNG setting will be
  used.}

  \item{rng.kind}{type of random number generator for
  bootstrap sampling. Only used when rng.seed is specified,
  see set.seed for details.}

  \item{iter.max}{maximum number of iterations for weighted
  Deming iterative algorithm.}

  \item{threshold}{numerical tolerance for weighted Deming
  iterative algorithm convergence.}

  \item{na.rm}{remove measurement pairs that contain
  missing values (Default is FALSE).}

  \item{NBins}{number of bins used when
  'reg.method="PaBaLarge"' to classify each slope in one of
  'NBins' bins covering the range of all slopes}
}
\value{
"MCResult" object containing regression results. The
function \code{\link{getCoefficients}} or
\code{\link{printSummary}} can be used to obtain or print a
summary of the results.  The function \code{\link{getData}}
allows to see the original data.  An S4 object of class
"MCResult" containing at least the following slots:
\item{data}{measurement data in wide format, one pair of
observations per sample. Includes samples ID, reference
measurement, test measurement.} \item{para}{numeric matrix
with estimates for slope and intercept, corresponding
standard deviations and confidence intervals.}
\item{mnames}{character vector of length two containing
names of analytical methods.} \item{regmeth}{type of
regression type used for parameter estimation.}
\item{cimeth}{method used for calculation of confidence
intervals.} \item{error.ratio}{ratio between squared
measurement errors of reference and test method, necessary
for Deming regression.} \item{alpha}{confidence level using
for calculation of confidence intervals.}
}
\description{
\code{mcreg} is used to compare two measurement methods by
means of regression analysis. Available methods comprise
ordinary and weighted linear regression, Deming and
weighted Deming regression and Passing-Bablok regression.
Point estimates of regression parameters are computed
together with their standard errors and confidence
intervals.
}
\details{
The regression analysis yields regression coefficients
'Inercept' and 'Slope' of the regression \eqn{Testmethod =
Intercept + Slope * Referencemethod}. There are methods for
computing the systematical bias between reference and test
method at a decision point Xc, \eqn{Bias(Xc) = Intercept +
(Slope-1) * Xc}, accompanied by its corresponding standard
error and confidence interval. One can use plotting method
\code{plotBias} for a comprehensive view of the
systematical bias.

Weighted regression for heteroscedastic data is available
for linear and Deming regression and implemented as a data
point weighting with the inverted squared value of the
reference method. Therefore calculation of weighted
regression (linear and Deming) is available only for
positive values (>0). Passing-Bablok regression is only
available for non-negative values (>=0).

Confidence intervals for regression parameters and bias
estimates are calculated either by using analytical methods
or by means of resampling methods ("jackknife",
"bootstrap", "nested bootstrap"). An analytical method is
available for all types of regression except for weighted
Deming. For Passing-Bablok regression the option
"analytical" calculates confidence intervals for the
regression parameters according to the non-parametric
approach given in the original reference. Note that
estimating resampling based confidence intervals for
Passing-Bablok regressions can take very long for larger
data sets due to the high computational complexity of the
algorithm. To mitigate this drawback an adaption of the
Passing-Bablok algorithm has been implemented
("PaBaLarge"), which yields approximative results. This
approach does not build the complete upper triangular
matrix of all 'n*(n-1)/2' slopes. It subdivides the range
of slopes into 'NBins' classes, and sorts each slope into
one of these bins. The remaining steps are the same as for
the exact "PaBa" algorithm, except that these are performed
on the binned slopes instead of operating on the matrix of
slopes.

The "jackknife" (or leave one out resampling) method was
suggested by Linnet for calculating confidence intervals of
regression parameters of Deming and weighted Deming
regression. It is possible to calculate jackknife
confidence intervals for all types of regression. Note that
we do not recommend this method for Passing-Bablok since it
has a tendency of underestimating the variability
(jackknife is known to yield incorrect estimates for errors
of quantiles).

The bootstrap method requires additionally choosing a value
for \code{method.bootstrap.ci}. If bootstrap is the method
of choice, "BCa", t-bootstrap ("tBoot") and simple
"quantile" confidence intervals are recommended (See Efron
B. and Tibshirani R.J.(1993),Carpenter J., Bithell J.
(2000)). The "nestedbootstrap" method can be very
time-consuming but is necessary for calculating t-bootstrap
confidence intervals for weighted Deming or Passing-Bablok
regression. For these regression methods there are no
analytical solutions for computing standard errors, which
therefore have to be obtained by nested bootstrapping.
}
\examples{
library("mcr")
    data(creatinine,package="mcr")
    x <- creatinine$serum.crea
    y <- creatinine$plasma.crea
    # Deming regression fit.
    # The confidence intercals for regression coefficients
    # are calculated with analytical method
    model1<- mcreg(x,y,error.ratio=1,method.reg="Deming", method.ci="analytical",
                     mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE)
    # Results
    printSummary(model1)
    getCoefficients(model1)
    plot(model1)
    # Deming regression fit.
    # The confidence intervals for regression coefficients
    # are calculated with bootstrap (BCa) method
    model2<- mcreg(x,y,error.ratio=1,method.reg="Deming",
                   method.ci="bootstrap", method.bootstrap.ci = "BCa",
                   mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE)
    compareFit(model1, model2)
}
\author{
Ekaterina Manuilova
\email{ekaterina.manuilova@contractors.roche.com}, Andre
Schuetzenmeister \email{andre.schuetzenmeister@roche.com},
Fabian Model \email{fabian.model@roche.com}
}
\references{
Bland, J. M., Altman, D. G. (1986) Statistical methods for
assessing agreement between two methods of clinical
measurement.  \emph{Lancet}, \bold{i:} 307--310.

Linnet, K. (1993) Evaluation of Regression Procedures for
Methods Comparison Studies.  \emph{CLIN. CHEM.}
\bold{39/3}, 424--432.

Linnet, K. (1990) Estimation of the Linear Relationship
between the Measurements of two Methods with Proportional
Errors.  \emph{STATISTICS IN MEDICINE}, Vol. \bold{9},
1463--1473.

Neter, J., Wassermann, W., Kunter, M. (1985) \emph{Applied
Statistical Models.} Richard D. Irwing, INC.

Looney, S. W. (2010) Statistical Methods for Assessing
Biomarkers.  \emph{Methods in Molecular Biology}, vol.
\bold{184}: \emph{Biostatistical Methods}. Human Press INC.

Passing, H., Bablok, W. (1983) A new biometrical procedure
for testing the equality of measurements from two different
analytical methods.  Application of linear regression
procedures for method comparison studies in clinical
chemistry, Part I.  \emph{J Clin Chem Clin Biochem}.  Nov;
\bold{21(11)}:709--20.

Efron, B., Tibshirani, R.J. (1993) \emph{An Introduction to
the Bootstrap}. Chapman and Hall.

Carpenter, J., Bithell, J. (2000) Bootstrap confidence
intervals: when, which, what? A practical guide for medical
statisticians.  \emph{Stat Med}, \bold{19 (9)}, 1141--1164.

\emph{CLSI EP9-A2}. Method Comparison and Bias Estimation
Using Patient Samples; Approved Guideline.
}
\seealso{
\code{\link{plotDifference}}, \code{\link{plot.mcr}},
\code{\link{getResiduals}}, \code{\link{plotResiduals}},
\code{\link{calcResponse}}, \code{\link{calcBias}},
\code{\link{plotBias}}, \code{\link{compareFit}}
}

