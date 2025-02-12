\name{show.particle.data}
\alias{show.particle.data}
\alias{show,particle.data-method}
\title{Show method for particle data class}
\description{
Display method for objects of class \code{particle.data}, which prints a summary of the key parameters and characteristics of the particle system based on its data type (simulation or experimental).
}
\usage{
\S4method{show}{particle.data}(object)
}
\arguments{
\item{object}{An object of class \code{particle.data}.}
}
\details{
This method displays essential information about the particle system. The output varies based on the data type:

For simulation data:
\itemize{
\item Type of particle interaction model
\item Number of time steps
\item Number of particles in the system
\item Noise variance parameter (sigma_0)
\item Interaction radius between particles
}

For experimental data:
\itemize{
\item Number of time steps
\item Average number of particles across all time steps
}
}
\references{
Vicsek, T., Czirok, A., Ben-Jacob, E., Cohen, I., & Shochet, O. (1995). \emph{Novel type of phase transition in a system of self-driven particles}, \emph{Physical Review Letters}, \bold{75}(6), 1226.

Chat'{e}, H., Ginelli, F., Gr'{e}goire, G., Peruani, F., & Raynaud, F. (2008). \emph{Modeling collective motion: variations on the Vicsek model}, \emph{The European Physical Journal B}, \bold{64}(3), 451-456.

Fang, X., & Gu, M. (2024). \emph{The inverse Kalman filter}. arXiv:2407.10089.
}
\seealso{
\code{\link{particle.data-class}} for the full class description,
\code{\link{simulate_particle}} for creating simulation data objects,
\code{\link{trajectory_data}} for creating experimental data objects
}
\keyword{methods}