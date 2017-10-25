\name{WaveletFitting}
\alias{WaveletFitting}
\title{Wavelet transform using Maximal overlap discrete wavelet transform (MODWT) algorithm}
\description{
Transforms the time series data by using hybrid MODWT algorithm using 'haar' filter.
}
\usage{
WaveletFitting(ts,Wvlevels,bndry,FFlag)
}
\arguments{
  \item{ts}{univariate time series}
  \item{Wvlevels}{The level of wavelet decomposition}
  \item{bndry}{The boundary condition of wavelet decomposition}
  \item{FFlag}{The FastFlag condition of wavelet decomposition: True or False}
     }
\value{
    \item{WaveletFitting}{The wavelet trasnform of the series}
  }

\references{
Percival D. B. and Walden A. T. 2000. Wavelet Methods for Time-Series Analysis. Cambridge Univ. Press, U.K.

Paul R. K., Prajneshu and Ghosh H. 2013. Wavelet Frequency Domain Approach for Modelling and Forecasting of Indian Monsoon Rainfall Time-Series Data. Journal of the Indian society of agricultural statistics, 67, 319 to 327.

Paul, R.K. and Birthal, P.S. 2015. Investigating rainfall trend over India using wavelet technique. Journal of Water and Climate Change, 7, 365 to 378.

Paul, R. K. 2015. ARIMAX-GARCH-WAVELET Model for forecasting volatile data. Model Assisted Statistics and Application, 10, 243 to252.
}
\examples{

N <- 100
PHI <- 0.2
THETA <- 0.1
SD <- 1
M <- 0
D <- 0.2
Seed <- 123

set.seed(Seed)
Sim.Series <- fracdiff::fracdiff.sim(n = N, ar = c(PHI), ma = c(THETA),
                                     d = D, rand.gen = rnorm, sd = SD, mu = M)
simts <- as.ts(Sim.Series$series)
Waveletlevels <- floor(log(length(simts))) # to obtain the maximum level for wavelet decomposition
WS <- WaveletFitting(ts=simts,Wvlevels=Waveletlevels,bndry='periodic',FFlag=TRUE)$WaveletSeries
}
