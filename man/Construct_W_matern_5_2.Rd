\name{Construct_W_matern_5_2}
\alias{Construct_W_matern_5_2}
\title{ 
The conditional covariance matrix for  matern covariance with roughness parameter 2.5
}
\description{
The conditional covariance matrix of the state in the dynamic linear model when kernel is the matern covariance with roughness parameter 2.5.
}
\usage{
Construct_W_matern_5_2(sigma2,delta_x,lambda,W0)
}
\arguments{
  \item{sigma2}{the variance parameter.
}
  \item{delta_x}{the distance between the sorted input.
}
  \item{lambda}{the transformed range parameter. }
    \item{W0}{the covariance matrix of the stationary distribution of the state. }

}


\value{
W matrix. 
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
