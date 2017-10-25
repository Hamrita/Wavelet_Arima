\name{WaveletFittingarma}
\alias{WaveletFittingarma}
\title{Wavelet-ARIMA hybrid model for forecasting}
\description{
Fits the time series data by using hybrid Wavelet-ARIMA algorithm.
}
\usage{
WaveletFittingarma(ts,Waveletlevels,boundary,FastFlag,MaxARParam,MaxMAParam,NForecast)
}
\arguments{
  \item{ts}{univariate time series}

  \item{Waveletlevels}{The level of wavelet decomposition}
  \item{boundary}{The boundary condition of wavelet decomposition}
  \item{FastFlag}{The FastFlag condition of wavelet decomposition: True or False}
  \item{MaxARParam}{The maximum AR order for auto.arima}
  \item{MaxMAParam}{The maximum MA order for auto.arima}
  \item{NForecast}{The forecast horizon: A positive integer}
   }
\value{
    \item{WaveletFittingarma}{The forecast of the series}
  }

\references{
Aminghafari, M. and Poggi, J.M. 2007. Forecasting time series using wavelets. Internationa Journal of Wavelets, Multiresolution and Inforamtion Processing, 5, 709 to 724

Percival D. B. and Walden A. T. 2000. Wavelet Methods for Time-Series Analysis. Cambridge Univ. Press, U.K.

Paul R. K., Prajneshu and Ghosh H. 2013. Wavelet Frequency Domain Approach for Modelling and Forecasting of Indian Monsoon Rainfall Time-Series Data. Journal of the Indian society of agricultural statistics, 67, 319 to 327.

Paul, R.K. and Birthal, P.S. 2015. Investigating rainfall trend over India using wavelet technique. Journal of Water and Climate Change, 7, 365 to 378.

Paul, R. K. 2015. ARIMAX-GARCH-WAVELET Model for forecasting volatile data. Model Assisted Statistics and Application, 10, 243 to 252.
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
#Waveletlevels <- floor(log(length(simts))) # to obtain the maximum level for wavelet decomposition
WaveletForecast<-WaveletFittingarma(ts=simts,Waveletlevels=floor(log(length(simts))),
boundary='periodic',FastFlag=TRUE,MaxARParam=5,MaxMAParam=5,NForecast=5)
}
