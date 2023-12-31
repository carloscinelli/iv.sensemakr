% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iv_functions.R
\name{iv_fit}
\alias{iv_fit}
\title{Instrumental Variable Estimation using the Anderson-Rubin approach}
\usage{
iv_fit(y, d, z, x = NULL, h0 = 0, alpha = 0.05)
}
\arguments{
\item{y}{\code{\link{numeric}} \code{\link{vector}} with the outcome.}

\item{d}{\code{\link{numeric}} \code{\link{vector}} with the treatment.}

\item{z}{\code{\link{numeric}} \code{\link{vector}} with the instrument.}

\item{x}{(optional) \code{\link{numeric}} \code{\link{matrix}} with observed covariates.}

\item{h0}{null hypothesis for the target parameter (the IV estimate).}

\item{alpha}{significance test for hypothesis tests and confidence intervals.}
}
\value{
An object of class \code{iv_fit}, containing:
\describe{
 \item{ \code{data} }{A \code{\link{data.frame}} with the data used for fitting the models.}
 \item{ \code{models} }{A \code{\link{list}} with the \code{\link{lm}} models used for obtaining the IV estimates. This includes the first-stage (FS), reduced-form (RF), and Anderson-Rubin (AR) regressions.}
 \item{ \code{estimates} }{A \code{\link{list}} with the summary information of IV estimates, as well as summary information of the auxiliary estimates of the FS, RF, and AR regression. }
 \item{\code{pars}}{A \code{\link{list}} with the parameters of the call, such as the null hypothesis \code{h0} and the significance level \code{alpha}.}
 }

 @returns
 The function returns an object of class \code{iv_fit}.
}
\description{
\code{iv_fit} computes instrumental variable estimates and confidence intervals using the Anderson-Rubin (AR) approach (Anderson and Rubin, 1949).  This approach is numerically identical to Fieller's theorem (Fieller, 1954). See Cinelli and Hazlett (2023) for further discussion.

The AR point estimate is numerically identical to the point estimate of two-stage least squares (2SLS) and it is given by the ratio of the reduced-form to the first-stage regression coefficient. Confidence intervals, however, are constructed differently. 2SLS is equivalent to using the delta-method to obtain the variance of the ratio estimator, and then proceeding by assuming the ratio is asymptotically normal. This approximation can fail when instruments are "weak." The Anderson-Rubin approach instead uses a test inversion procedure to construct confidence intervals. This procedure has correct coverage regardless of instrument strength, at the (inevitable) cost of eventually obtaining unbounded confidence intervals.
}
\references{
Anderson, T.W. and Rubin, H. (1949), Estimation of the parameters of a single equation in a complete system of stochastic equations, Annals of Mathematical Statistics, 20, 46-63.

Fieller, E. C. (1954). Some problems in interval estimation. Journal of the Royal Statistical Society: Series B (Methodological), 16(2), 175-185.

Cinelli, C. and Hazlett, C. (2023), "An Omitted Variable Bias Framework for Sensitivity Analysis of Instrumental Variables."
}
