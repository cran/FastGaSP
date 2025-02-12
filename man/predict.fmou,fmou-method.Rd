\name{predict.fmou}
\alias{predict.fmou}
\alias{predict.fmou,fmou-method}
\title{ 
Prediction and uncertainty quantification on the future observations using a FMOU model.}


\description{
This function predicts the future observations using a FMOU model. Uncertainty quantification is available. 
}


\usage{
\S4method{predict.fmou}{fmou}(object, param_est, step=1, interval=FALSE, interval_data=TRUE)
}


\arguments{
  \item{object}{an objecft of  class \code{fmou.}}
  \item{param_est}{estimated parameters in the FMOU model. Obtained by the results of \code{fit.fmou()}. }
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

Lin, Y., Liu, X., Segall, P., & Gu, M. (2025). Fast data inversion for high-dimensional dynamical systems from noisy measurements. arXiv preprint arXiv:2501.01324.

}

\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}


\examples{

## generate simulated data
library(FastGaSP)
library(rstiefel)

d = 5  # number of latent factors
k = 20 # length of observation at each time step
n = 500 # number time step
noise_level = 1 # variance of noise

U = rustiefel(k, k) # factor loading matrix
z = matrix(NA, d, n)
sigma_2 = runif(d, 0.5, 1)
rho = runif(d, 0.95, 1)
for(l in 1:d){
  R = matrix(NA, n, n)
  diag(R) = 1
  for(ir in 1:n){
    for(ic in 1:n){
      R[ir, ic] = rho[l]^(abs(ir-ic)) * R[ir, ir]
    }
  }
  R = (sigma_2[l]/(1-rho[l]^2) )* R
  z[l, ] = t(chol(R)) \%*\% rnorm(n)
}

signal = U[,1:d] \%*\% z
y = signal + matrix(rnorm(n*k,mean=0,sd=sqrt(noise_level)),k,n)

##constucting the fmou.model
fmou.model=fmou(output=y, d=d, est_U0=TRUE, est_sigma0_2=TRUE)

## estimate the parameters
em_alg <- fit.fmou(fmou.model, M=500)

## two-step-ahead prediction
pred_2step <- predict.fmou(fmou.model,em_alg, step=c(1:2))

}
