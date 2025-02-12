\name{predict.gppca}
\alias{predict.gppca}
\alias{predict.gppca,gppca-method}
\title{ 
Prediction and uncertainty quantification on the future observations using GPPCA.}

\description{
This function predicts the future observations using a GPPCA model. Uncertainty quantification is available. 
}

\usage{
\S4method{predict.gppca}{gppca}(object, param, A_hat, step=1, interval=FALSE, interval_data=TRUE)
}

\arguments{
  \item{object}{an object of  class \code{gppca.}}
  \item{param}{estimated parameters (beta, sigma2, sigma0_2).}
  \item{A_hat}{estimated factor loading matrix.}
  \item{step}{a vector. Number of steps to be predicted. Default is 1. } 
  \item{interval}{a bool value, default is \code{FALSE}. If \code{TRUE}, the 95\% predictive intervals are computed.}
  \item{interval_data}{a bool value, default is \code{TRUE}. If \code{TRUE}, the 95\% predictive intervals of the observations are computed. Otherwise, the 95\% predictive intervals of the mean of the observation are computed.} 
}

\value{
  \item{pred_mean}{the predictive mean. }
  \item{pred_interval_95lb}{the 95\% lower bound of the interval.}
  \item{pred_interval_95ub}{the 95\% upper bound of the interval.}
}

\references{

Gu, M., & Shen, W. (2020), Generalized probabilistic principal component analysis of correlated data, \emph{Journal of Machine Learning Research, 21}(13), 1-41.

}

\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}


\examples{

library(FastGaSP)
library(rstiefel)

matern_5_2_funct <- function(d, beta_i) {
  cnst <- sqrt(5.0)
  matOnes <- matrix(1, nrow = nrow(d), ncol = ncol(d))
  result <- cnst * beta_i * d
  res <- (matOnes + result + (result^2) / 3) * exp(-result)
  return(res)
}

n=200
k=8
d=4

beta_real=0.01
sigma_real=1
sigma_0_real=sqrt(.01)
input=seq(1,n,1)
R0_00=as.matrix(abs(outer(input,input,'-')))
R_r = matern_5_2_funct(R0_00, beta_real)
L_sample = t(chol(R_r))
kernel_type='matern_5_2'


input=sort(input)
delta_x=input[2:length(input)]-input[1:(length(input)-1)]
A=rustiefel(k, d)  ##sample from Stiefel manifold
Factor=matrix(0,d,n)
for(i in 1: d){
  Factor[i,]=sigma_real^2*L_sample\%*\%rnorm(n)
}
output=A\%*\%Factor+matrix(rnorm(n*k,mean=0,sd=sigma_0_real),k,n)
  
##constucting the gppca.model
gppca_obj <- gppca(input, output, d, shared_params = TRUE,est_d=FALSE)
## estimate the parameters
gppca_fit <- fit.gppca(gppca_obj)

## two-step-ahead prediction
param <- c(gppca_fit$est_beta, gppca_fit$est_sigma2, gppca_fit$est_sigma0_2)
gppca_pred <- predict.gppca(gppca_obj, param, gppca_fit$est_A, step=1:3)
gppca_pred$pred_mean

}
