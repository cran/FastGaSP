\name{Get_R_y}
\alias{Get_R_y}
\title{ 
The multiplication of R with y
}
\description{
This function computes the product of the R matrix and the output vector, where R is the correlation matrix for a dynamic linear model (DLM). Instead of explicitly forming the Cholesky decomposition of R, this function computes the product as \eqn{L (L^T y)}, where \eqn{L} is the Cholesky decomposition of R. This is achieved using the forward filtering algorithm for efficient computation.
}
\usage{
Get_R_y(GG, Q, K, output)
}
\arguments{
  \item{GG}{a list of matrices defined in the dynamic linear model.}
  
  \item{Q}{a vector defined in the dynamic linear model.} % Each element is the variance of the corresponding state.
  
  \item{K}{a matrix defined in the filtering algorithm for the dynamic linear model.}
  
  \item{output}{a vector of observations.}
}
\value{
A vector representing the product of the R matrix and the output vector, where \eqn{R} is the correlation matrix for a dynamic linear model.
}
\references{
Hartikainen, J. and Sarkka, S. (2010). \emph{Kalman filtering and smoothing solutions to temporal Gaussian process regression models}. \emph{Machine Learning for Signal Processing (MLSP), 2010 IEEE International Workshop}, 379-384.

Fang, X., & Gu, M. (2024). \emph{The inverse Kalman filter}. arXiv:2407.10089.

M. Gu, Y. Xu (2019), \emph{Fast nonseparable Gaussian stochastic process with application to methylation level interpolation}. \emph{Journal of Computational and Graphical Statistics}, In Press, arXiv:1711.11501.

Campagnoli P, Petris G, Petrone S. (2009), \emph{Dynamic linear model with R}. Springer-Verlag New York.

}
\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}
\seealso{
\code{\link{Get_Q_K}} for more details on \eqn{K} and \eqn{Q} matrices, 
\code{\link{Get_L_inv_y}} for \eqn{L^{-1}y} computation, 
\code{\link{Get_L_t_y}} for \eqn{L^T y} computation, 
\code{\link{Get_L_y}} for \eqn{L y} computation, 
\code{\link{Get_L_t_inv_y}} for \eqn{(L^T)^{-1}y} computation.
}
\keyword{internal}
