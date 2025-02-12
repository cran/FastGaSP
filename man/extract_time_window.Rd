\name{extract_time_window}
\alias{extract_time_window}
\title{Extract time window from particle data}
\description{
Extracts a specified time window from a particle.data object while preserving all relevant tracking information and parameters. Works with both simulation and experimental data.
}
\usage{
extract_time_window(data_obj, first_frame, last_frame)
}
\arguments{
\item{data_obj}{An object of class \code{\linkS4class{particle.data}}.}
\item{first_frame}{Integer specifying the first frame to include (must be >= 1).}
\item{last_frame}{Integer specifying the last frame to include (must be less than the total number of frames).}
}
\value{
Returns a new \code{\linkS4class{particle.data}} object containing:
\describe{
\item{px_list, py_list}{Position data for the selected time window}
\item{vx_list, vy_list}{Velocity data for the selected time window}
\item{theta_list}{Angle data if present in original object}
\item{particle_tracking}{Tracking information for the selected frames (experimental data)}
\item{data_type}{Original data type ("simulation" or "experiment")}
\item{n_particles}{Particle counts (constant for simulation, time series for experimental)}
\item{T_time}{Number of time steps in the extracted window}
\item{model, sigma_0, radius}{Original model parameters (for simulation data)}
\item{D_y}{Original dimension of the output space}
}
}