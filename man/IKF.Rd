\name{IKF}
\alias{IKF}
\title{
Inverse Kalman Filter - The multiplication of R with y with given kernel type
}
\description{
This function computes the product of the R matrix and the output vector for a given kernel type, where R is the correlation matrix for a dynamic linear model (DLM). Instead of explicitly forming the Cholesky decomposition of R, this function computes the product as \eqn{L (L^T y)}, where \eqn{L} is the Cholesky decomposition of R. This is achieved using the forward filtering algorithm for efficient computation.
}
\usage{
IKF(beta, tilde_nu, delta_x, output, kernel_type='matern_5_2')
}

\arguments{
  \item{beta}{A scalar representing the inverse range parameter in the kernel function.}
  \item{tilde_nu}{A numerical value representing the nugget or the stabilizing parameter.}
  \item{delta_x}{A numeric vector of successive differences of inputs.}
  \item{output}{The output vector for which the product with R is computed.}
  \item{kernel_type}{A string specifying the kernel type. Must be either "matern_5_2" or "exp".}
}

\value{
A vector representing the product of the R matrix and the output vector
}
\references{
Fang, X., & Gu, M. (2024). \emph{The inverse Kalman filter}. arXiv:2407.10089.
}
\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}

\seealso{
\code{\link{Get_R_y}} for computing the product of R and y with G, Q, and K matrices.
}
\examples{
# Helper function for Matern 2.5 correlation
matern25_correlation <- function(d, beta) {
  #' Compute Matern 2.5 correlation matrix given distance matrix and beta parameter
  #'
  #' @param d Distance matrix or array of any dimensions
  #' @param beta Inverse range parameter 
  #'
  #' @return Correlation matrix with same dimensions as d
  #'
  #' @details
  #' Formula: k(d) = (1 + sqrt(5)*beta*d + 5*beta^2*d^2/3) * exp(-sqrt(5)*beta*d)
  
  # Compute sqrt(5)*beta*d
  sqrt5_beta_d <- sqrt(5) * beta * d
  
  # Compute the correlation
  R <- (1 + sqrt5_beta_d + (5 * beta^2 * d^2) / 3) * exp(-sqrt5_beta_d)
  
  return(R)
}

# Example: Comparing IKF with direct matrix computation
n <- 2000                           # number of observations
input <- sort(runif(n))             # sorted input points
delta_x <- input[-1] - input[-n]    # consecutive differences
u <- rnorm(n)                       # random output vector
beta <- 10                          # inverse range parameter

# Test 1: Noise-free scenario
# Non-robust IKF (no stabilization parameter)
x_non_robust_IKF <- IKF(beta = beta, tilde_nu = 0, delta_x = delta_x, 
                        output = u, kernel_type = 'matern_5_2')

# Robust IKF (with stabilization parameter)
tilde_nu <- 0.1  # stabilization parameter
x_IKF <- IKF(beta = beta, tilde_nu = tilde_nu, delta_x = delta_x, 
             output = u, kernel_type = 'matern_5_2') - tilde_nu * u

# Direct matrix computation for comparison
R0 <- abs(outer(input, input, '-'))  # distance matrix
R <- matern25_correlation(R0, beta)  # correlation matrix
x_direct <- R \%*\% u          

# Compare results (should be nearly identical)
cat("Maximum difference (non-robust IKF vs direct):", 
    max(abs(x_direct - x_non_robust_IKF)), "\n")
cat("Maximum difference (robust IKF vs direct):", 
    max(abs(x_direct - x_IKF)), "\n")

# Test 2: With noise
V <- 0.2  # noise variance
x_IKF_noisy <- IKF(beta = beta, tilde_nu = V, delta_x = delta_x, 
                   output = u, kernel_type = 'matern_5_2')

# Direct computation with noise
x_direct_noisy <- (R + V * diag(n)) \%*\% u

# Compare results
cat("Maximum difference (IKF vs direct):", 
    max(abs(x_direct_noisy - x_IKF_noisy)), "\n")

}
