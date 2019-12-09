\name{coef.lod_lm}
\alias{coef.lod_lm}
\title{
  Extract lod_lm Coefficients
}
\description{
  Extracts estimates regression coefficients from object of class "\code{lod_lm}".
}
\usage{
  ## S3 method for class 'lod_lm'
  \method{coef}{lod_lm}(object, ...)
}
\arguments{

  \item{object}{An object of class "\code{lod_lm}", usually, a result of a call to \code{lod_lm}}

  \item{...}{further arguments passed to or from other methods.}

  }

\value{
  Coefficients extracted from \code{object} as a named numeric vector.
}
\author{
Kevin Donovan, kmdono02@ad.unc.edu.

Maintainer: Kevin Donovan <kmdono02@ad.unc.edu>
}
\references{
  \insertRef{may}{lodr}
}

\seealso{
  \code{\link[lodr]{fitted.lod_lm}} and \code{\link[lodr]{residuals.lod_lm}} for related methods; \code{\link[lodr]{lod_lm}} for model fitting.

  The generic functions \code{\link[stats]{fitted}} and \code{\link[stats]{residuals}}.
}
\examples{
 library(lodr)
 ## Using example dataset provided in lodr package: lod_data_ex
 ## 3 covariates: X1, X2, X3 with X2 and X3 subject to a lower limit of
 ## detection of 0

 fit <- lod_lm(data=lod_data_ex, frmla=Y~X1+X2+X3, lod=c(0,0),
                   var_LOD=c("X2", "X3"))
 coef(fit)
}