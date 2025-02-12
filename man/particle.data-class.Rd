\name{particle.data-class}
\alias{particle.data-class}
\alias{particle.data}
\docType{class}
\title{Particle trajectory data class}
\description{
S4 class for storing and analyzing particle trajectory data from both simulations and experimental observations. This class supports different models including Vicsek and can handle both position and velocity data along with optional angle information and particle tracking capabilities for experimental data.
}
\section{Objects from the Class}{
Objects of this class can be created in two ways:
\itemize{
\item For simulation data: Using \code{\link{simulate_particle}} that computes particle trajectories under physical models
\item For experimental data: Using \code{\link{trajectory_data}} to save particle trajectories while handling varying numbers of particles between time steps
}
}
\section{Slots}{
    \describe{
        \item{\code{px_list}:}{Object of class \code{list}. List of x-positions at each time step.}
        \item{\code{py_list}:}{Object of class \code{list}. List of y-positions at each time step.}
        \item{\code{vx_list}:}{Object of class \code{list}. List of x-velocities at each time step.}
        \item{\code{vy_list}:}{Object of class \code{list}. List of y-velocities at each time step.}
        \item{\code{theta_list}:}{Object of class \code{listOrNULL}. Optional list of particle velocity angles at each time step.}
        \item{\code{particle_tracking}:}{Object of class \code{listOrNULL}. List of data frames containing particle mappings between consecutive frames (primarily for experimental data).}
        \item{\code{data_type}:}{Object of class \code{character}. Type of data: either "simulation" or "experiment".}
        \item{\code{n_particles}:}{Object of class \code{numeric}. Number of particles (constant for simulation data, or a vector recording the number of particles at each time step for experimental data).}
        \item{\code{T_time}:}{Object of class \code{numeric}. Total number of time steps.}
        \item{\code{D_y}:}{Object of class \code{numeric}. Dimension of the output space.}
        \item{\code{model}:}{Object of class \code{characterOrNULL}. Type of particle interaction model (e.g., "Vicsek"). NULL for experimental data.}
        \item{\code{sigma_0}:}{Object of class \code{numericOrNULL}. Noise variance parameter used in the model. NULL for experimental data.}
        \item{\code{radius}:}{Object of class \code{numericOrNULL}. Interaction radius between particles. NULL for experimental data.}
    }
}

\section{Methods}{
\describe{
  \item{show:}{Method for displaying summary information about the particle.data object.}
  \item{fit:}{Method for fitting the latent factor model to data using the IKF-CG algorithm, which returns a \code{particle.est} object containing estimated parameters and predictions. See \code{\link{fit.particle.data}} for detailed documentation.}
}
}



\references{
Vicsek, T., Czirok, A., Ben-Jacob, E., Cohen, I., & Shochet, O. (1995). \emph{Novel type of phase transition in a system of self-driven particles}, \emph{Physical Review Letters}, \bold{75}(6), 1226.

Chat'{e}, H., Ginelli, F., Gr'{e}goire, G., Peruani, F., & Raynaud, F. (2008). \emph{Modeling collective motion: variations on the Vicsek model}, \emph{The European Physical Journal B}, \bold{64}(3), 451-456.
}

\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}

\seealso{
\code{\link{simulate_particle}} for simulating particle trajectories,
\code{\link{trajectory_data}} for saving experimantal particle trajectories,
\code{\link{fit.particle.data}} for model fitting methods
}
\keyword{classes}