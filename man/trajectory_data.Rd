\name{trajectory_data}
\alias{trajectory_data}
\title{Convert experimental particle tracking data to particle.data object}
\description{
Processes experimental particle tracking data and creates a standardized particle.data object. This function handles time series data with varying numbers of particles across time steps and maintains particle identity tracking between consecutive frames.
}
\usage{
trajectory_data(particle_data)
}
\arguments{
\item{particle_data}{A data frame containing particle tracking data with the following required columns:
\itemize{
\item \code{px}, \code{py}: Particle positions in x and y coordinates
\item \code{vx}, \code{vy}: Particle velocities in x and y directions
\item \code{time}: Time step identifier (integer). Must be consecutive integers starting from the minimum time value
\item \code{particleID}: Unique particle identifier for tracking across frames
}
}
}
\value{
Returns an S4 object of class \code{\linkS4class{particle.data}} with:
  \describe{
    \item{px_list, py_list}{Lists of particle x and y positions at each time step}
    \item{vx_list, vy_list}{Lists of particle x and y velocities at each time step}
    \item{theta_list}{List of particle angles computed from velocities}
    \item{particle_tracking}{List of data frames containing particle mappings between consecutive frames}
    \item{data_type}{"experiment"}
    \item{n_particles}{Vector recording number of particles at each time step}
    \item{T_time}{Total number of time steps}
    \item{D_y}{Dimension of the output space (set to 1)}
  }
}

\examples{
# Create sample tracking data
sample_data <- data.frame(time = rep(1:3, each = 5),
                          particleID = rep(1:5, 3),
                          px = rnorm(15),py = rnorm(15),
                          vx = rnorm(15),vy = rnorm(15)
)
# Convert to particle.data object
traj <- trajectory_data(sample_data)
# Display summary
show(traj)
}
\seealso{
\code{\link{particle.data-class}} for details on the particle.data class structure
}


