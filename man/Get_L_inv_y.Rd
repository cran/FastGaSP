\name{Get_L_inv_y}
\alias{Get_L_inv_y}
\title{ 
The multiplication of the inverse of L with y
}
\description{
This function computes the product of the inverse of the L matrix and the output vector, where the L matrix is the Cholesky decomposition of the correlation matrix. Instead of computing the Cholesky matrix, we compute it using the forward filtering algorithm. 
}
\usage{
Get_L_inv_y(GG,Q,K,output)
}
\arguments{
  \item{GG}{a list of matrices defined in the dynamic linear model. 
}
%%     ~~Describe \code{response} here~~
    \item{Q}{a vector defined in the dynamic linear model.    }
    \item{K}{a matrix defined in the filtering algorithm for the dynamic linear model.  }
    \item{output}{a vector of output.}

}


\value{
A vector representing the product of the inverse of the L matrix and the output vector, where the L matrix is the Cholesky decomposition of the correlation matrix. 
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
\seealso{\code{\link{Get_Q_K}} for more details about \eqn{Q} vector and \eqn{K} matrix,
\code{\link{Get_L_t_y}} for \eqn{L^T y} computation, 
\code{\link{Get_L_y}} for \eqn{L y} computation, 
\code{\link{Get_L_t_inv_y}} for \eqn{(L^T)^{-1}y} computation.
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
