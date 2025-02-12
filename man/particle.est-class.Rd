\name{particle.est-class}
\alias{particle.est-class}
\alias{particle.est}
\docType{class}
\title{Particle interaction estimation class}
\description{S4 class for storing estimated parameters and predictions for particle interaction models. }

\section{Objects from the Class}{
Objects of this class are created by the \code{\link{fit.particle.data}} (via \code{fit}) method when applied to \code{\linkS4class{particle.data}} objects to estimate interaction parameters and make predictions.
}

\section{Slots}{
    \describe{
        \item{\code{data_type}:}{Object of class \code{character}. Specifies the type of data ("simulation" or "experiment").}
        \item{\code{model}:}{Object of class \code{characterOrNULL}. Specifies the model type for simulation data (e.g., "Vicsek" or "two_interactions_Vicsek"). NULL for experimental data.}
        \item{\code{D_y}:}{Object of class \code{numeric}. Dimension of the output space.}
        \item{\code{num_interaction}:}{Object of class \code{numeric}. Number of interactions.}
        \item{\code{parameters}:}{Object of class \code{numeric}. Vector of estimated parameters with length 2*D_y + 1:
          \itemize{
          \item First D_y elements: beta (inverse range parameters)
          \item Next D_y elements: tau (variance-noise ratios)
          \item Last element: interaction radius
          }
        }
        \item{\code{sigma_2_0_est}:}{Object of class \code{numeric}. Estimated noise variance.}
        \item{\code{predictions}:}{object of class \code{listOrNULL}. Contains predicted means and 95\% confidence intervals (lower and upper bounds) for the particle interactions if testing inputs are given.}
        \item{\code{training_data}:}{Object of class \code{list}. Contains the training data used in the GP model, obtained using the estimated interaction radius.}
        \item{\code{gp_weights}:}{Object of class \code{matrix}. Contains the weights from the GP computation (A^T_j Sigma_y^(-1) y) used for prediction, with each column corresponding to a type of interaction j.}
    }
}

\section{Methods}{
    \describe{
        \item{show:}{Method for displaying summary information about the estimated parameters.}
    }
}

\references{
Fang, X., & Gu, M. (2024). \emph{The inverse Kalman filter}. arXiv:2407.10089.
}

\seealso{
  \code{\link{fit.particle.data}} for more details about how to create a \code{particle.est} object.
  \code{\link{particle.data-class}} for the input data structure
}

\keyword{classes}