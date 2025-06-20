\name{A_times_x_particle}
\alias{A_times_x_particle}
\title{Matrix-vector multiplication for particle systems}
\description{
  Performs the matrix-vector multiplication A*x for particle systems, where A is a sparse matrix stored only the non-zero entries
}

\usage{
A_times_x_particle(output, A_all_v, num_neighbors_vec, D, N)
}

\arguments{
  \item{output}{A numeric vector containing the input vector for multiplication.}
  
  \item{A_all_v}{A numeric vector containing the interaction matrices in vectorized form.}
  
  \item{num_neighbors_vec}{An integer vector specifying the number of neighbors for each particle.}
  
  \item{D}{An integer specifying the dimension of the output vector per particle.}
  
  \item{N}{An integer specifying the total dimension of the output vector.}
}


\value{
  Returns a numeric vector containing the result of the matrix-vector multiplication.
}

\references{
Fang, X., & Gu, M. (2024). \emph{The inverse Kalman filter}. arXiv:2407.10089.
}

\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}

\keyword{internal}