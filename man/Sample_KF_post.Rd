\name{Sample_KF_post}
\alias{Sample_KF_post}
\title{ 
Sample the posterior distribution of the process using the backward smoothing algorithm  
}
\description{
This function samples the posterior distribution of the process using the backward smoothing algorithm. 
}
\usage{
Sample_KF_post(index_obs, C_R_K_Q,W0,GG,W,VV,output,kernel_type,sample_type)
}
\arguments{
    \item{index_obs}{a vector  where the entries with 1 have observations and entries with 0 have no observation.}
        \item{C_R_K_Q}{a list of matrices to compute the inverse covariance matrix in the dynamic linear model.}

  \item{GG}{a list of matrices defined in the dynamic linear model. 
}
  \item{W}{a list of coefficient matrices defined in the dynamic linear model. 
}
  \item{C0}{the covariance matrix of the stationary distribution defined in the dynamic linear model. 
}

  \item{VV}{a numerical value of the variance of the nugget parameter. }
  \item{output}{a vector of the output. }

  \item{kernel_type}{a \code{character} to specify the type of kernel to use. The current version supports kernel_type to be "matern_5_2" or "exp", meaning that the matern kernel with roughness parameter being 2.5 or 0.5 (exponent kernel), respectively. }
    \item{sample_type}{a integer to specify the type of sample we need. 0 means the states. 1 means the first value of each state vector. 2 means the noisy observations. }

}


\value{
A matrix of the posterior samples. 
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
\seealso{\code{\link{Sample_KF_post}} for more details about sampling from the posterior distribution.
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
