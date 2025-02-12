\name{gppca-class}
\docType{class}
\alias{gppca-class}
\title{GPPCA class}
\description{
An S4 class for generalized probabilistic principal component analysis of correlated data.
}

\section{Objects from the Class}{
Objects of this class are created and initialized using the \code{\link{gppca}} function to set up the estimation.
}

\section{Slots}{
	 \describe{
	  \item{\code{input}:}{object of class \code{vector}, the length is equivalent to the number of observations.}
    \item{\code{output}:}{object of class \code{matrix}. The observation matrix.}
    \item{\code{d}:}{object of class \code{integer} to specify the number of latent factors.}
    \item{\code{est_d}:}{object of class \code{logical}, default is \code{FALSE}. If \code{TRUE}, d will be estimated by either variance matching (when noise level is given) or information criteria (when noise level is unknown). Otherwise, d is fixed, and users must assign a value to \code{d}.}
    \item{\code{shared_params}:}{object of class \code{logical}, default is \code{TRUE}. If \code{TRUE}, the latent processes share the correlation and variance parameters. Otherwise, each latent process has distinct parameters.}
    \item{\code{kernel_type}:}{a \code{character} to specify the type of kernel to use. The current version supports kernel_type to be "matern_5_2" or "exponential", meaning that the matern kernel with roughness parameter being 2.5 or 0.5 (exponent kernel), respectively. }
    }
}

\section{Methods}{
  \describe{
    \item{fit.gppca}{See \code{\link[=fit.gppca]{fit.gppca}} for details.}
    \item{predict.gppca}{See \code{\link[=predict.gppca]{predict.gppca}} for details.}
	 }
}

\references{

Gu, M., & Shen, W. (2020), Generalized probabilistic principal component analysis of correlated data, \emph{Journal of Machine Learning Research, 21}(13), 1-41.

}


\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}

\seealso{ \code{\link{gppca}} for more details about how to create a \code{gppca} object.}


