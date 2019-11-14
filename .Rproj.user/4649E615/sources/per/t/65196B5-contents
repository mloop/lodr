\name{lod_lm}
\alias{lod_lm}
\title{
  Fitting Linear Models with Covariates Subject to a Limit of Detection
}
\description{
  \code{lod_lm} is used to fit linear models while taking into account limits of detection for corresponding covariates.  It carries out the method detailed in \insertCite{may;textual}{lodr} with regression coefficient standard errors calculated using bootstrap resampling.
}
\usage{
  lod_lm(data, frmla, lod, var_LOD,
  nSamples = 250, convergenceCriterion = 0.001, boots = 25)
}
\arguments{
  \item{data}{a required data frame (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables in the model.  If not specified, a corresponding error is returned.}

  \item{frmla}{an object of class "\link[stats]{formula}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.}

  \item{lod}{a numeric vector (or object coercible by \link[base]{as.numeric}) specifying the limit of detection for each covariates specified in \code{var_LOD} (in the same order as the covariates in \code{var_LOD}).}

  \item{var_LOD}{a character vector specifying which covariates in the model (\code{frmla}) are subject to limits of detection.}

  \item{nSamples}{an integer specifying the number of samples to generate for each subject with covariate values outside of their limits of detection.  For more details, see \insertCite{may;textual}{lodr}.  The default is 250.}

  \item{convergenceCriterion}{a number specifying the smallest difference between iterations required for the regression coefficient estimation process to complete.  The default is 0.001.}

  \item{boots }{a number specifying the number of bootstrap resamples used for the standard error estimation process for the regression coefficient estimates.  The default is 25.}

  }
\details{
  Models for \code{lod_lm} are specified the same as models for \code{lm}.  A typical model as the form \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a series of terms separated by \code{+} which specifies a linear predictor for \code{response}.  A formula has an implied intercept term.

  In the dataset used with lod_lm, values outside of the limits of detection need to be denoted by the value of the lower limit of detection.  Observations with values marked as missing by \code{NA} are removed by the model fit procedure as done with \code{lm}.
}
\value{
  \code{lod_lm} returns an object of \link[base]{class}) "\code{lod_lm}".  The function \code{summary} prints a summary of the results in the same format as with an object of \link[base]{class}) "\code{lm}".  The generic accessor functions \code{coef}, \code{fitted} and \code{residuals} extract various useful features of the value returned by \code{lod_lm}.

  An object of \link[base]{class}) "\code{lod_lm}" is a list containing the following components:

  \item{coefficients}{a named vector of regression coefficient estimates.}
  \item{boot_SE}{a named vector of regression coefficient estimate bootstrap standard error estimates.}
  \item{fitted.values}{the fitted mean values for subjects with covariates within their limits of detection.}
  \item{rank}{the numeric rank of the fitted linear model}
  \item{residuals}{the residuals, that is response minus fitted values, for subjects with covariates within their limits of detection.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{model}{the model frame used.}
  \item{call}{the matched call.}
  \item{terms}{the \code{\link[stats]{terms}} object used.}
}
\author{
Kevin Donovan, kmdono02@ad.unc.edu.

Maintainer: Kevin Donovan <kmdono02@ad.unc.edu>
}
\references{
  \insertRef{may}{lodr}
}

\seealso{
  \code{\link[lodr]{summary.lod_lm}} for summaries of the results from \code{lod_lm}

  The generic functions \code{link[stats]{coef}}, \code{link[stats]{fitted.values}} and \code{link[stats]{residuals}}.
}
\examples{
     library(lodr)
     ## Using example dataset provided in lodr package: lod_data_ex
     ## 3 covariates: X1, X2, X3 with X2 and X3 subject to a lower limit of
     ## detection of 0

     fit <- lod_lm(data=lod_data_ex, frmla=Y~X1+X2+X3, lod=c(0,0),
                   var_LOD=c("X2", "X3"))
     summary(fit)
}
