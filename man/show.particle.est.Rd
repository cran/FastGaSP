\name{show.particle.est}
\alias{show,particle.est-method}
\title{Show method for particle estimation class}
\description{
  Display method for objects of class \code{particle.est}, which prints a summary of the estimated parameters for the particle interaction model.
}
\usage{
  \S4method{show}{particle.est}(object)
}
\arguments{
  \item{object}{An object of class \code{particle.est}.}
}
\details{
This method displays essential information about the estimated model parameters, including:
  \itemize{
    \item Data type (simulation or experiment)
    \item Model type (for simulation data only)
    \item Dimension of output space
    \item Estimated parameters:
      \itemize{
      \item beta parameters (inverse range)
      \item tau parameters (variance-noise ratio)
      \item interaction radius
    }
  }
}
\references{
Fang, X., & Gu, M. (2024). \emph{The inverse Kalman filter}. arXiv:2407.10089.
}

\seealso{
\code{\link{particle.est-class}} for details of the particle estimation class
}
