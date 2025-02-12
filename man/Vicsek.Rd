\name{Vicsek}
\alias{Vicsek}
\title{Vicsek Model Simulation}
\description{
  Simulates particle movement according to the Vicsek model, where particles align their velocities with neighboring particles within a specified radius, subject to noise. The model implements both Gaussian and uniform noise options.
}
\usage{
Vicsek(p0, v0, theta0, v_abs, n_t, T_sim, h, cut_r, 
       sigma_0, noise_type = "Gaussian")
}
\arguments{
  \item{p0}{A numeric vector of initial positions for all particles, structured as (x1, y1, x2, y2, ..., xn, yn).}
  \item{v0}{A numeric vector of initial velocities for all particles, structured similarly to p0.}
  \item{theta0}{A numeric vector of initial angles for all particles.}
  \item{v_abs}{A numeric value specifying the absolute velocity (speed) of particles.}
  \item{n_t}{An integer specifying the number of particles.}
  \item{T_sim}{An integer specifying the number of time steps to simulate.}
  \item{h}{A numeric value specifying the time step size.}
  \item{cut_r}{A numeric value specifying the interaction radius within which particles align.}
  \item{sigma_0}{A numeric value specifying the noise strength.}
  \item{noise_type}{A character string specifying the type of noise: either "Gaussian" (default) or "Uniform".}
}
\value{
  Returns a list with three components:
  \item{pos}{A matrix of dimension (2*n_t) x (T_sim+1) containing particle positions at each time step.}
  \item{v}{A matrix of dimension (2*n_t) x (T_sim+1) containing particle velocities at each time step.}
  \item{theta}{A matrix of dimension n_t x (T_sim+1) containing particle angles at each time step.}
}
\references{
Vicsek, T., Czirok, A., Ben-Jacob, E., Cohen, I., & Shochet, O. (1995). \emph{Novel type of phase transition in a system of self-driven particles}, \emph{Physical Review Letters}, \bold{75}(6), 1226.
}

\keyword{internal}