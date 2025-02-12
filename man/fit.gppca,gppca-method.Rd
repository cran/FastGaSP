\name{fit.gppca}
\alias{fit.gppca}
\alias{fit.gppca,gppca-method}
\title{ 
Parameter estimation for generalized probabilistic principal component analysis of correlated data.}

\description{
This function estimates the parameters for generalized probabilistic principal component analysis of correlated data.
}

\usage{
\S4method{fit.gppca}{gppca}(object, sigma0_2=NULL, d_ub=NULL)
}

\arguments{
  \item{object}{an object of  class \code{gppca.}}
  \item{sigma0_2}{variance of noise. Default is \code{NULL}. User should specify a value when it is known in real data.}
  \item{d_ub}{upper bound of d when d is estimated. Default is \code{NULL}.}
}

\value{
  \item{est_A}{the estimated factor loading matrix.}
  \item{est_beta}{the estimated inversed range parameter.}
  \item{est_sigma0_2}{the estimated variance of noise.}
  \item{est_sigma2}{the estimated variance parameter.}
  \item{mean_obs}{the predictive mean of the observations.}
  \item{mean_obs_95lb}{the lower bound of the 95\% posterior credible intervals of predictive mean.}
  \item{mean_obs_95ub}{the upper bound of the 95\% posterior credible intervals of predictive mean.}
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

## MSE between predictive mean of observations and true mean
Y_mean=A\%*\%Factor
mean((gppca_fit$mean_obs-Y_mean)^2)
sd(Y_mean)

## coverage of 95\% posterior credible intervals
sum(gppca_fit$mean_obs_95lb<=Y_mean & gppca_fit$mean_obs_95ub>=Y_mean)/(n*k)

## length of 95\% posterior credible intervals
mean(abs(gppca_fit$mean_obs_95ub-gppca_fit$mean_obs_95lb))

}
