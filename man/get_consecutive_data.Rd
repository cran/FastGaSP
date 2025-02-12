\name{get_consecutive_data}
\alias{get_consecutive_data}
\title{Extract consecutive time steps data from particle trajectories}
\description{
  Extracts paired data from consecutive time steps for a specified variable in particle trajectory data. This function is particularly useful for experimental data, where it uses particle tracking information to maintain particle identity between frames.
}
\usage{
get_consecutive_data(data_obj, variable = c("vx", "vy", "px", "py", "theta"))
}
\arguments{
  \item{data_obj}{An object of class \code{\linkS4class{particle.data}} containing particle trajectory data.}
  \item{variable}{A character string specifying which variable to extract. Must be one of "vx" (x-velocity), "vy" (y-velocity), "px" (x-position), "py" (y-position), or "theta" (angle).}
}
\value{
  Returns a list with two components:
  \item{start}{A list of length T_time containing the data at each time step t.}
  \item{end}{A list of length T_time containing the corresponding data at time step t+1.}
  
  For each time t, start[[t]] and end[[t]] contain paired measurements for the same particles at consecutive time steps.
}
\examples{
\dontrun{
vx_pairs = get_consecutive_data(data_obj, "vx")
vx_list = vx_pairs$start
vx_end_list = vx_pairs$end
}
}

\keyword{internal}