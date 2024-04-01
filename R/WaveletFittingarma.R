#'@title Wavelet Transform Using Maximal Overlap Discrete Wavelet Transform (MODWT) Algorithm
#' @description Transforms the time series data by using hybrid MODWT algorithm.
#' @param ts Univariate time series
#' @param Wvlevels The level of wavelet decomposition
#' @param WFilter  Wavelet filter use in the decomposition
#' @param bndry The boundary condition of wavelet decomposition:'periodic' or 'reflection'
#' @param FFlag The FastFlag condition of wavelet decomposition: True or False
#' @import fracdiff forecast stats wavelets
#'
#' @return
#' \itemize{
#'   \item WaveletSeries - The wavelet trasnform of the series
#' }
#' @export
#'
#' @examples
#' data<-rnorm(100,mean=100,sd=50)
#' WaveletFitting(ts=data,Wvlevels=3,WFilter='haar',bndry='periodic',FFlag=TRUE)
#' @references
#' \itemize{
#'\item Aminghafari, M. and Poggi, J.M. 2007. Forecasting time series using wavelets. Internationa Journal of Wavelets, Multiresolution and Inforamtion Processing, 5, 709 to 724

#' \item Percival D. B. and Walden A. T. 2000. Wavelet Methods for Time-Series Analysis. Cambridge Univ. Press, U.K.

#' \item Paul R. K., Prajneshu and Ghosh H. 2013. Wavelet Frequency Domain Approach for Modelling and Forecasting of Indian Monsoon Rainfall Time-Series Data. Journal of the Indian society of agricultural statistics, 67, 319 to 327.
#' }
WaveletFitting <- function(ts,WFilter="haar",Wvlevels,bndry='periodic',FFlag=TRUE)
{
  mraout <- wavelets::mra(ts, filter=WFilter, n.levels=Wvlevels,boundary=bndry, fast=FFlag, method="modwt")
  WaveletSeries <- cbind(do.call(cbind,mraout@D),mraout@S[[Wvlevels]])
  return(list(WaveletSeries=WaveletSeries,WVSeries=mraout))
}

#'@title Wavelet-ARIMA hybrid model for forecasting
#' @description Fits the time series data by using hybrid Wavelet-ARIMA algorithm.
#' @param ts univariate time series
#' @param filter Wavelet filter use in the decomposition
#' @param Waveletlevels The level of wavelet decomposition
#' @param boundary The boundary condition of wavelet decomposition
#' @param FastFlag The FastFlag condition of wavelet decomposition: True or False
#' @param MaxARParam The maximum AR order for auto.arima
#' @param MaxMAParam The maximum MA order for auto.arima
#' @param NForecast The forecast horizon: A positive integer
#' @import fracdiff forecast stats wavelets
#' @return
#' \itemize{
#'   \item Finalforecast - Forecasted value
#'   \item FinalPrediction - Predicted value of train data
#' }
#' @export
#'
#' @examples
#' N <- 100
#'PHI <- 0.2
#'THETA <- 0.1
#' SD <- 1
#'M <- 0
#'D <- 0.2
#' Seed <- 123

#'set.seed(Seed)
#'Sim.Series <- fracdiff::fracdiff.sim(n = N,ar=c(PHI),ma=c(THETA),d=D,rand.gen =rnorm,sd=SD,mu=M)
#'simts <- as.ts(Sim.Series$series)
#'WaveletForecast<-WaveletFittingarma(ts=simts,filter ='la8',Waveletlevels=floor(log(length(simts))),
#'MaxARParam=5,MaxMAParam=5,NForecast=5)
#' @references
#' \itemize{
#'\item Aminghafari, M. and Poggi, J.M. 2012. Nonstationary time series forecasting using wavelets and kernel smoothing. Communications in Statistics-Theory and Methods, 41(3),485-499.
#' \item Paul, R.K. A and Anjoy, P. 2018. Modeling fractionally integrated maximum temperature series in India in presence of structural break. Theory and Applied Climatology 134, 241–249.
#' }
WaveletFittingarma<- function(ts,filter="haar",Waveletlevels,boundary='periodic',FastFlag=TRUE, MaxARParam,MaxMAParam,NForecast)

{
  WS <- WaveletFitting(ts=ts,WFilter=filter,Wvlevels=Waveletlevels,bndry=boundary,FFlag=FastFlag)$WaveletSeries
  AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL
  #-----------------------------------------------------------#
  # Fitting of ARIMA model to the Wavelet Coef                #
  #-----------------------------------------------------------#
  for(WVLevel in 1:ncol(WS))
  {
    ts <- NULL
    ts <- WS[,WVLevel]
    WaveletARMAFit <- forecast::auto.arima(x=as.ts(ts), d=NA, D=NA, max.p=MaxARParam, max.q=MaxMAParam,stationary=FALSE,
                                           seasonal=FALSE,ic=c("aic"), allowdrift=FALSE, allowmean=TRUE,stepwise = TRUE)
    WaveletARIMAPredict <- WaveletARMAFit$fitted
    WaveletARIMAForecast <- forecast::forecast(WaveletARMAFit,h=NForecast)
    AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletARIMAPredict)
    AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletARIMAForecast$mean))
  }
  Finalforecast <- rowSums(AllWaveletForecast,na.rm = T)
  FinalPrediction <- rowSums(AllWaveletPrediction,na.rm = T)
  return(list(Finalforecast=Finalforecast,FinalPrediction=FinalPrediction))
}
