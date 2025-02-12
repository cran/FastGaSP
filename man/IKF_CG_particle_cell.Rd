\name{IKF_CG_particle_cell}
\alias{IKF_CG_particle_cell}
\title{Inverse Kalman Filter with Conjugate Gradient for Particle Systems}
\description{
  Implements the IKF-CG (Inverse Kalman Filter - Conjugate Gradient) algorithm for the one-interaction physical model with 1D output, different particle numbers for different time, and non-identity diagonal noise. 
}
\usage{
IKF_CG_particle_cell(param, kernel_type, delta_x_all, output, A_all_v, 
                     sort_d_all_ix, sigma_2_vec, num_neighbors_vec, 
                     tilde_nu, D, n_t_record, tol = 1e-6, maxIte = 1000)
}
\arguments{
  \item{param}{A numeric vector containing model parameters. The first element is the log of beta (inverse range parameter), and the second element is the log of tau (variance ratio parameter).}
  \item{kernel_type}{A string specifying the covariance kernel type: "matern_5_2" or "exp".}
  \item{delta_x_all}{A numeric vector of successive differences of sorted inputs.}
  \item{output}{A numeric vector of observations.}
  \item{A_all_v}{A numeric vector containing the non-zero entries in matrix A.}
  \item{sort_d_all_ix}{An integer vector of indices for sorting distances.}
  \item{sigma_2_vec}{A numeric vector of variances for each time point.}
  \item{num_neighbors_vec}{An integer vector specifying number of neighbors for each observation.}
  \item{tilde_nu}{A numeric value for numerical stabilization.}
  \item{D}{An integer specifying the dimension of observations.}
  \item{n_t_record}{An integer vector containing number of particles at each time point.}
  \item{tol}{A numeric value specifying convergence tolerance (default: 1e-6).}
  \item{maxIte}{An integer specifying maximum number of iterations (default: 1000).}
}
\value{
  Returns a list containing:
  \item{x}{A numeric vector of solved coefficients.}
  \item{resid}{A numeric vector of residuals at each iteration.}
  \item{ite}{An integer indicating the number of iterations performed.}
}
\references{
  Fang, X., & Gu, M. (2024). \emph{The inverse Kalman filter}. arXiv:2407.10089.
}

\keyword{internal}