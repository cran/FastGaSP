\name{fit.fmou}
\alias{fit.fmou}
\alias{fit.fmou,fmou-method}
\title{ 
The fast EM algorithm of multivariate Ornstein-Uhlenbeck processes}


\description{
This function implements an efficient EM algorithm to estimate the parameters in the FMOU model, a latent factor model with a fixed or estimated orthogonal factor loading matrix, where each latent factor is modeled as an O-U (Ornstein-Uhlenbeck) process.
}


\usage{
\S4method{fit.fmou}{fmou}(object, M=50, threshold=1e-4,
         track_iterations=FALSE,track_neg_log_lik=FALSE,
         U_init=NULL, rho_init=NULL, sigma2_init=NULL, d_ub=NULL)
}


\arguments{
  \item{object}{an objecft of  class \code{fmou.}}
  \item{M}{number of iterations in the EM algorithm, default is 50.}
  \item{threshold}{stopping criteria with respect to predictive mean of observations, default is 1e-4.}
  \item{track_iterations}{a bool value, default is \code{FALSE}. If \code{TRUE}, the estimations in each EM iteration will be recorded and returned. }
  \item{track_neg_log_lik}{a bool value, default is \code{FALSE}. If \code{TRUE}, the negative log marginal likelihood of the output in each EM iteration will be recorded and returned.}
  \item{U_init}{user-specified initial factor loading matrix in the EM algorithm. Default is \code{NULL}. The dimension is k*d, where k is the length of observations at each time step and d is the number of latent factors.}
  \item{rho_init}{user-specified initial correlation parameters in the EM algorithm. Default is \code{NULL}. The length is equal to the number of latent factors.}
  \item{sigma2_init}{user-specified initial variance parameters in the EM algorithm. Default is \code{NULL}. The length is equal to the number of latent factors.}
  \item{d_ub}{upper bound of d when d is estimated. Default is null.}
}

\value{
  \item{output}{the observation matrix.}
  \item{U}{the estimated (or fixed) factor loading matrix.}
  \item{post_z_mean}{the posterior mean of latent factors.}
  \item{post_z_var}{the posterior variance of latent factors.}
  \item{post_z_cov}{the posterior covariance between two consecutive time steps of a latent process.}
  \item{mean_obs}{the predictive mean of the observations.}
  \item{mean_obs_95lb}{the lower bound of the 95\% posterior credible intervals of predictive mean.}
  \item{mean_obs_95ub}{the upper bound of the 95\% posterior credible intervals of predictive mean.}
  \item{sigma0_2}{estimated variance of noise.}
  \item{rho}{estimated correlation parameters.}
  \item{sigma2}{estimated variance parameters}
  \item{num_iterations}{number of iterations in the EM algorithm.}
  \item{d}{the estimated (or fixed) number of latent factors.}
  \item{record_sigma0_2}{estimated variance of noise in each iteration.}
  \item{record_rho}{estimated correlation parameters in each iteration.}
  \item{record_sigma2}{estimation variance parameters in each iteration.}
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

## root mean square error (RMSE) of predictive mean of observations
sqrt(mean((em_alg$mean_obs-signal)^2))

## standard deviation of (truth) mean of observations
sd(signal)

## estimated variance of noise
em_alg$sigma0_2

}
