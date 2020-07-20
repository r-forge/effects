\name{EffectMethods}
\alias{effSources}
\alias{effSources.default}
\alias{effSources.gls}
\alias{effSources.glmmPQL}
\alias{effSources.rlmerMod}
\alias{effSources.clm}
\alias{effSources.clm2}
\alias{effSources.clmm}
\alias{effSources.betareg}


\title{Functions For Constructing Effect Displays for Many Modeling Paradigms}
\description{The \code{effects} package can be used with many different modelling types beyond the basic models described at \code{\link{Effect}}.  The \code{effSources} function is used to tell the \code{Effect} method where to find the information to compute effects.  Details are given in the vignette \href{../doc/methods-supported-by-effects.pdf}{Regression Models Supported by the Effects Package}.  Unless you are trying to extend \code{effects} to a different regression model, you will have no use for this function.  
}

\usage{

\method{effSources}{default}(mod)

\method{effSources}{gls}(mod)

\method{effSources}{glmmPQL}(mod)

\method{effSources}{rlmerMod}(mod)

\method{effSources}{clm}(mod)

\method{effSources}{clm2}(mod)

\method{effSources}{clmm}(mod)

\method{effSources}{betareg}(mod)

}

\arguments{
  \item{mod}{a fitted model object of the appropriate class. The default method works for many regression models for which there is no specific method.}
  }

\details{
The \code{Effect} method, and its equivalents, expect a regression object to provide certain quantities like model formula, coefficient estimates and the variance matrix of the estimates in standard places. The \code{effSources} methods return a list of needed quantities if they are not in standard places.}

\value{See the vignette \href{../doc/methods-supported-by-effects.pdf}{Regression Models Supported by the Effects Package}.}




\author{John Fox \email{jfox@mcmaster.ca}, Sanford Weisberg \email{sandy@umn.edu}}


%\examples{
%}

\keyword{internal}
