\name{Get_Q_K}
\alias{Get_Q_K}
\title{ 
 one-step-ahead predictive variance and Kalman gain
}
\description{
This function computes the one-step-ahead predictive variance and Kalman gain.
}
\usage{
Get_Q_K(GG,W,C0,VV)
}
\arguments{
%%     ~~Describe \code{response} here~~
    \item{GG}{a list of matrices defined in the dynamic linear model.}
    \item{W}{a list of matrices defined in the dynamic linear model.}
    \item{C0}{a matrix defined in the dynamic linear model.}
    \item{VV}{a numerical value for the nugget.}
}


\value{
A list of 2 items for Q and K. 
}

\references{
%% ~put references to the literature/web site here ~

Hartikainen, J. and Sarkka, S. (2010). \emph{Kalman filtering and smoothing solutions to temporal gaussian process regression models}.  \emph{Machine Learning for Signal Processing (MLSP), 2010 IEEE International Workshop},  379-384.

M. Gu, Y. Xu (2019), \emph{fast nonseparable gaussian stochastic process with application to methylation level interpolation}.  \emph{Journal of Computational and Graphical Statistics}, In Press, arXiv:1711.11501.

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
