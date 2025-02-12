\name{IKF_CG_particle}
\alias{IKF_CG_particle}
\title{IKF-CG algorithm for one-interaction physical model with 1D output}
\description{
  Implements the IKF-CG (Inverse Kalman Filter - Conjugate Gradient) algorithm for the one-interaction physical model with 1D output. This function provides an efficient computational method for large-scale particle systems.
}

\usage{
IKF_CG_particle(param, kernel_type, delta_x_all, output, A_all_v, 
                sort_d_all_ix, num_neighbors_vec, tilde_nu,
                D, N, tol = 1e-6, maxIte = 1000)
}

\arguments{
  \item{param}{A numeric vector containing model parameters. The first element is the log of beta (inverse range parameter), and the second element is the log of tau (variance ratio parameter).}
  
  \item{kernel_type}{A string specifying the kernel type. Must be either "matern_5_2" or "exp".}
  
  \item{delta_x_all}{A numeric vector of successive differences of sorted inputs.}
  
  \item{output}{A numeric vector representing the output values.}
  
  \item{A_all_v}{A numeric vector containing the non-zero entries in matrix A.}
  
  \item{sort_d_all_ix}{An integer vector containing sorted indices for distances.}
  
  \item{num_neighbors_vec}{An integer vector specifying the number of neighbors for each particle.}
  
  \item{tilde_nu}{A numeric value representing the stabilizing parameter.}
  
  \item{D}{An integer specifying the dimension of the output vector per particle.}
  
  \item{N}{An integer specifying the total dimension of the output vector.}
  
  \item{tol}{A numeric value specifying the convergence tolerance. Default is 1e-6.}
  
  \item{maxIte}{An integer specifying the maximum number of iterations. Default is 1000.}
}


\value{
  Returns a list containing three elements:
  \item{x}{A numeric vector containing the estimated state variables.}
  \item{resid}{A numeric vector containing the residuals at each iteration.}
  \item{ite}{An integer indicating the number of iterations performed.}
}



\references{
Fang, X., & Gu, M. (2024). \emph{The inverse Kalman filter}. arXiv:2407.10089.
}



\keyword{internal}