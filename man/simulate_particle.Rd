\name{simulate_particle}
\alias{simulate_particle}
\title{Simulate particle trajectories}
\description{
  Simulates particle trajectories using either the standard Vicsek model or a two-interaction variation of the Vicsek model.
}
\usage{
  simulate_particle(v_abs, n_t = 100, T_sim = 5, h = 0.1,
    cut_r = 0.5, sigma_0 = 0.1,
    noise_type = "Gaussian", model = "Vicsek")
}
\arguments{
  \item{v_abs}{Absolute velocity magnitude for all particles.}
  \item{n_t}{Number of particles (default: 100).}
  \item{T_sim}{Total simulation time steps (default: 5).}
  \item{h}{Time step size for numerical integration (default: 0.1).}
  \item{cut_r}{Radius of interaction between particles (default: 0.5).}
  \item{sigma_0}{Standard deviation of noise (default: 0.1).}
  \item{noise_type}{Distribution of noise: "Gaussian" or "Uniform" (default: "Gaussian").}
  \item{model}{Type of interaction model: "Vicsek", "unnormalized_Vicsek", or "two_interactions_Vicsek" (default: "Vicsek").}
}
\value{
Returns an S4 object of class \code{\linkS4class{particle.data}}. See \code{\link{particle.data-class}} for details of the returned object structure.
}

\references{
Vicsek, T., Czirok, A., Ben-Jacob, E., Cohen, I., & Shochet, O. (1995). \emph{Novel type of phase transition in a system of self-driven particles}, \emph{Physical Review Letters}, \bold{75}(6), 1226.

Chat'{e}, H., Ginelli, F., Gr'{e}goire, G., Peruani, F., & Raynaud, F. (2008). \emph{Modeling collective motion: variations on the Vicsek model}, \emph{The European Physical Journal B}, \bold{64}(3), 451-456.

Fang, X., & Gu, M. (2024). \emph{The inverse Kalman filter}. arXiv:2407.10089.
}

\examples{
#--------------------------------------------------
# Example: Simulate using standard Vicsek model
#--------------------------------------------------
vx_abs=0.5
vy_abs=0.5
v_abs=sqrt(vx_abs^2+vy_abs^2)
sim1 <- simulate_particle(v_abs=v_abs)

#--------------------------------------------------
# Example: Simulate using unnormalized variation
#--------------------------------------------------
sim2 <- simulate_particle(v_abs=v_abs, model = 'unnormalized_Vicsek')

#--------------------------------------------------
# Example: Simulate using two-interaction variation
#--------------------------------------------------
sim3 <- simulate_particle(v_abs=v_abs, model = 'two_interactions_Vicsek')
}
\seealso{
\code{\link{particle.data-class}} for details on the particle.data class structure
}
