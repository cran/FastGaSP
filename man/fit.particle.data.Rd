\name{fit.particle.data}
\alias{fit.particle.data}
\alias{fit,particle.data-method}
\title{Fit method for particle data}
\description{
Estimates interaction parameters for particle systems using trajectory data with the IKF-CG (Inverse Kalman Filter - Conjugate Gradient) approach. Supports both simulation and experimental data.
}
\usage{
  \S4method{fit}{particle.data}(
    object, param, cut_r_max=1, est_param = TRUE, nx=NULL, ny=NULL,
    kernel_type = "matern_5_2", tilde_nu = 0.1, tol = 1e-6,
    maxIte = 1000, output = NULL, ho_output = NULL, 
    testing_inputs=NULL, compute_CI = TRUE, num_interaction = (length(param)-1)/2,
    data_type = object@data_type, model = object@model, 
    apolar_vicsek = FALSE, direction = NULL
  )
}
\arguments{
  \item{object}{An object of class \code{particle.data} containing the trajectory data.}
  \item{param}{Numeric vector of parameters. Should contain 2*num_interaction + 1 elements: first num_interaction elements are log of inverse range parameters (beta), next num_interaction elements are log of variance-noise ratios (tau), and the final element is log(radius/(cut_r_max-radius)) where radius is the interaction radius.}
  \item{cut_r_max}{Numeric value specifying the maximum interaction radius to consider during estimation (default: 1).}
  \item{est_param}{If TRUE, param is used as initial values for parameter optimization. If FALSE, param is treated as fixed parameters for prediction (default: TRUE).}
  \item{nx}{An integer specifying the number of grid points along the x-axis (horizontal direction). If NULL, automatically calculated as floor((px_max-px_min)/cut_r_max), where px_max and px_min represent the maximum and minimum x-coordinates of all particles.}  
  \item{ny}{An integer specifying the number of grid points along the y-axis (vertical direction). If NULL, automatically calculated as floor((py_max-py_min)/cut_r_max), where py_max and py_min represent the maximum and minimum y-coordinates of all particles.}
  \item{kernel_type}{Character string specifying the kernel type: 'matern_5_2' (default) or 'exp'.}
  \item{tilde_nu}{Numeric value for stabilizing the IKF computation (default: 0.1).}
  \item{tol}{Numeric value specifying convergence tolerance for the conjugate gradient algorithm (default: 1e-6).}
  \item{maxIte}{Integer specifying maximum iterations for the conjugate gradient algorithm (default: 1000).}
  \item{output}{Numerical vector (default = NULL). Used for residual bootstrap when different outputs but same inputs are needed.}
  \item{ho_output}{Numerical vector (default = NULL). Used for residual bootstrap when different hold-out outputs but same inputs are needed.}
  \item{testing_inputs}{Matrix of inputs for prediction (NULL if only performing parameter estimation). Each column represents testing inputs for one interaction.}
  \item{compute_CI}{When TRUE, computes the 95\% credible interval for testing_inputs (default: TRUE).}
  \item{num_interaction}{Integer specifying number of interactions to predict (default: (length(param_ini)-1)/2).}
  \item{data_type}{Character string indicating data type ("simulation" or "experiment").}
  \item{model}{Character string specifying the model type (e.g., "unnormalized_Vicsek").}
  \item{apolar_vicsek}{When TRUE, considers only neighboring particles moving in the same direction (default: FALSE).}
  \item{direction}{Modeling direction ('x' or 'y') for experimental data analysis.}
}
\value{
Returns an object of class \code{\linkS4class{particle.est}}. See \code{\link{particle.est-class}} for details.
}
\references{
Fang, X., & Gu, M. (2024). \emph{The inverse Kalman filter}. arXiv:2407.10089.
}
\examples{
# Simulate data
vx_abs <- 0.5
vy_abs <- 0.5
v_abs <- sqrt(vx_abs^2+vy_abs^2)
sim <- simulate_particle(v_abs=v_abs, model = 'unnormalized_Vicsek')
show(sim)

# Set up testing inputs and initial parameters
testing_n <- 200
testing_inputs <- matrix(as.numeric(seq(-1,1,length.out=testing_n)),nr=1)
cut_r_max=1.5
param_ini <- log(c(0.3,1000,0.3/(cut_r_max-0.3)))  # Initial parameter values

# Fit model to simulation data
est <- fit(sim,param=param_ini,cut_r_max=1.5, testing_inputs = testing_inputs)
show(est)
}
\keyword{methods}