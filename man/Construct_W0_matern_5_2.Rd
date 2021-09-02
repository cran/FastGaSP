\name{Construct_W0_matern_5_2}
\alias{Construct_W0_matern_5_2}
\title{ 
covariance of the stationary distribution of the state when kernel is the Matern covariance with roughness parameter 2.5. 
}
\description{
This function computes covariance of the stationary distribution of the state when kernel is the Matern covariance with roughness parameter 2.5. 
}
\usage{
Construct_W0_matern_5_2(sigma2,lambda)
}
\arguments{
  \item{sigma2}{the variance parameter.
}
  \item{lambda}{the transformed range parameter. }
}


\value{
W0 matrix. 
}

\references{
%% ~put references to the literature/web site here ~

Hartikainen, J. and Sarkka, S. (2010). \emph{Kalman filtering and smoothing solutions to temporal gaussian process regression models}.  \emph{Machine Learning for Signal Processing (MLSP), 2010 IEEE International Workshop},  379-384.

M. Gu, Y. Xu (2019), \emph{fast nonseparable Gaussian stochastic process with application to methylation level interpolation}.  \emph{Journal of Computational and Graphical Statistics}, In Press, arXiv:1711.11501.

Campagnoli P, Petris G, Petrone S. (2009), \emph{Dynamic linear model with R}. Springer-Verlag New York.

}

\author{
\packageAuthor{FastGaSP}

Maintainer: \packageMaintainer{FastGaSP}
}

\keyword{internal}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
