\name{fmou-class}
\docType{class}
\alias{fmou-class}
\title{ FMOU class}
\description{
An S4 class for fast parameter estimation in the FMOU model, a latent factor model with a fixed or estimated orthogonal factor loading matrix, where each latent factor is modeled as an O-U (Ornstein-Uhlenbeck) process.
}

\section{Objects from the Class}{
Objects of this class are created and initialized using the \code{\link{fmou}} function to set up the estimation.
}

\section{Slots}{
	 \describe{
    \item{\code{output}:}{object of class \code{matrix}. The observation matrix.}
    \item{\code{d}}{object of class \code{integer} to specify the number of latent factors.}
    \item{\code{est_d}}{object of class \code{logical}, default is \code{FALSE}. If \code{TRUE}, d will be estimated by either variance matching (when noise level is given) or information criteria (when noise level is unknown). Otherwise, d is fixed, and users must assign a value to \code{d}.}
    \item{\code{est_U0}}{object of class \code{logical}, default is \code{TRUE}. If \code{TRUE}, the factor loading matrix (U0) will be estimated. Otherwise, U0 is fixed.}
    \item{\code{est_sigma0_2}}{object of class \code{logical}, default is \code{TRUE} . If \code{TRUE}, the variance of the noise will be estimated. Otherwise, it is fixed.}
    \item{\code{U0}}{object of class \code{matrix}. The fixed factor loading matrix. Users should assign a k*d matrix to it when \code{est_U0=False}. Here k is the length of observations at each time step.}
    \item{\code{sigma0_2}}{object of class \code{numeric}. Variance of noise. User should assign a value to it when \code{est_sigma0_2=False}.}
    }
}

\section{Methods}{
  \describe{
    \item{fit.fmou}{See \code{\link[=fit.fmou]{fit.fmou}}.}
    \item{predict.fmou}{See \code{\link[=predict.fmou]{predict.fmou}}.}
	 }
}

\references{

Lin, Y., Liu, X., Segall, P., & Gu, M. (2025). Fast data inversion for high-dimensional dynamical systems from noisy measurements. arXiv preprint arXiv:2501.01324.

}

\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}

\seealso{ \code{\link{fmou}} for more details about how to create a \code{fmou} object.}


