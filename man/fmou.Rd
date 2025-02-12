\name{fmou}
\alias{fmou}
\alias{fmou-method}
\title{ 
Setting up the FMOU model}
\description{
Creating an \code{fmou} class for fmou, a latent factor model with a fixed or estimated orthogonal factor loading matrix, where each latent factor is modeled as an O-U (Ornstein-Uhlenbeck) process.
}

\usage{
  fmou(output, d, est_d=FALSE, est_U0=TRUE, est_sigma0_2=TRUE, U0=NULL, sigma0_2=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{output}{a k*n observation matrix, where k is the length of observations at each time step and n is the number of time steps.}
    \item{d}{number of latent factors.}
    \item{est_d}{a bool value, default is \code{FALSE}. If \code{TRUE}, d will be estimated by either variance matching (when noise level is given) or information criteria (when noise level is unknown). Otherwise, d is fixed, and users must assign a value to argument d.}
    \item{est_U0}{a bool value, default is \code{TRUE}. If \code{TRUE}, the factor loading matrix (U0) will be estimated. Otherwise, U0 is fixed.}
    \item{est_sigma0_2}{a bool value, default is \code{TRUE} . If \code{TRUE}, the variance of the noise will be estimated. Otherwise, it is fixed.}
    \item{U0}{the fixed factor loading matrix. Users should assign a k*d matrix to it when \code{est_U0=FALSE}. }
    \item{sigma0_2}{variance of noise. User should assign a value to it when \code{est_sigma0_2=FALSE}.}

}

\value{
\code{fmou} returns an S4 object of class \code{fmou}.
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
n = 100 # number time step
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
}

