\name{f_Vicsek_variation}
\alias{f_Vicsek_variation}
\title{Modified Vicsek Interaction Function}
\description{
  A modified interaction function for the Vicsek model that defines the interaction strength between particles based on their distance. This function is used as one of the interactions in the two-interaction Vicsek model.
}
\usage{
f_Vicsek_variation(r, a = 0.02, b = 1, r_min = 0.01, r_max = 0.8, beta = 20)
}
\arguments{
  \item{r}{numeric, the distance between particles}
  \item{a}{numeric, strength parameter for the short-range interaction term. Default is 0.02}
  \item{b}{numeric, strength parameter for the linear term. Default is 1}
  \item{r_min}{numeric, minimum distance parameter for the interaction term. Default is 0.01}
  \item{r_max}{numeric, maximum distance parameter. Default is 0.8}
  \item{beta}{numeric, scaling parameter for the overall interaction. Default is 20}
}
\details{
  The function implements a modified Vicsek interaction with three components:
  \itemize{
    \item A short-range term (-a/(r+r_min))
    \item A linear term (-b*(r-r_max))
    \item A constant offset term (a/r_max)
  }
  The final value is scaled by the beta parameter.
}
\value{
  Returns a numeric value representing the interaction strength at distance r.
}
\examples{
# Calculate interaction values at various distances
r_seq <- seq(0.01, 1, by = 0.01)
interaction <- f_Vicsek_variation(r_seq)
plot(r_seq, interaction, type = "l", 
     xlab = "Distance", ylab = "Interaction Strength",
     main = "Vicsek Variation Interaction Function")
}

\references{
Chat'{e}, H., Ginelli, F., Gr'{e}goire, G., Peruani, F., & Raynaud, F. (2008). \emph{Modeling collective motion: variations on the Vicsek model}, \emph{The European Physical Journal B}, \bold{64}(3), 451-456.

Fang, X., & Gu, M. (2024). \emph{The inverse Kalman filter}. arXiv:2407.10089.

}

\keyword{internal}
