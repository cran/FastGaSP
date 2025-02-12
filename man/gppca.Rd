\name{gppca}
\alias{gppca}
\alias{gppca-method}
\title{Setting up the GPPCA model}
\description{Creating an \code{gppca} class for generalized probabilistic 
principal component analysis of correlated data.}

\usage{
gppca(input,output, d, est_d=FALSE, shared_params=TRUE, kernel_type="matern_5_2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{input}{a vector for the sorted input locations. The length is equal to the number of observations.}
    \item{output}{a k*d matrix for the observations at the sorted input locations. Here k is the number of locations and n is the number of observations.}
    \item{d}{number of latent factors.}
    \item{est_d}{a bool value, default is \code{FALSE}. If \code{TRUE}, d will be estimated by either variance matching (when noise level is given) or information criteria (when noise level is unknown). Otherwise, d is fixed, and users must assign a value to argument d.}
    \item{shared_params}{a bool value, default is \code{TRUE}. If \code{TRUE}, the latent processes share the correlation and variance parameters. Otherwise, each latent process has distinct parameters.}
    \item{kernel_type}{a \code{character} to specify the type of kernel to use. The current version supports kernel_type to be "matern_5_2" or "exponential", meaning that the matern kernel with roughness parameter being 2.5 or 0.5 (exponent kernel), respectively. }
    % description of inputs of the gppca() function
}

\value{
\code{gppca} returns an S4 object of class \code{gppca}.
}

\references{

Gu, M., & Shen, W. (2020), Generalized probabilistic principal component analysis of correlated data, \emph{Journal of Machine Learning Research, 21(13)}, 1-41.

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
output=A%*%Factor+matrix(rnorm(n*k,mean=0,sd=sigma_0_real),k,n)
  
##constucting the gppca.model
gppca_obj <- gppca(input, output, d, shared_params = TRUE,est_d=FALSE)

}

