\name{predictorEffects}
\alias{predictorEffect}
\alias{predictorEffect.svyglm}
\alias{predictorEffect.default}
\alias{predictorEffects}
\alias{predictorEffects.default}
\title{
Functions For Computing Predictor Effects
}
\description{
Alternatives to the \code{Effect} and \code{allEffects} functions that use a different paradigm for conditioning in an effect display.  The user specifies one predictor, either continuous or a factor, for the horizontal axis of a plot, and the function determines the appropriate plot to display (which is drawn by \code{plot}).
}
\usage{
predictorEffect(predictor, mod, xlevels=c(50, 5), ...)

\method{predictorEffect}{svyglm}(predictor, mod, xlevels=c(50, 5), ...)

\method{predictorEffect}{default}(predictor, mod, xlevels=c(50, 5), ...)

predictorEffects(mod, predictors, ...)

\method{predictorEffects}{default}(mod, predictors = ~ ., ...)

}
\arguments{
  \item{mod}{
A model object.  Supported models include all those described on the help page for \code{\link{Effect}}.
}
  \item{predictor}{quoted name of the focal predictor.}
  \item{predictors}{
If the default \code{~ .}, a predictor effect plot is drawn for each predictor (not regressor) in a model.  Otherwise, this should be a one-sided formula listing the first-order predictors for which predictor effect plots should be drawn.
}
  \item{xlevels}{this argument is used to set the number of levels (i.e., values) for any focal predictor that is not a factor. It is similar to the \code{xlevels} argument of \code{\link{Effect}}, except that it may be a numeric vector with two elements, the first of which applies to the focal numeric predictor and the second to
  any conditioning numeric predictors. The default is \code{c(50, 5)}, that is, 50 levels of the focal predictor and 5 levels of conditioning predictors. If a single number is given, it applies \emph{both} to the focal numeric predictor and conditioning numeric predictors. The
  \code{xlevels} argument may also be a list with named elements, as explained in the help for \code{\link{Effect}}.}
  \item{\dots}{
Additional arguments passed to \code{\link{Effect}}.
}
}
\details{
Effect plots view a fitted regression function E(Y|X) in (sequences of) two-dimensional plots using conditioning and slicing.  The functions described here use a different method of determining the conditioning and slicing than \code{allEffects} uses.  The predictor effect of a focal predictor, say \code{x1},is the  usual effect for the generalized interaction of \code{x1} with all the other predictors in a model.  When a predictor effect object is plotted, the focal predictor is by default plotted on the horizontal axis.

For example, in the model \code{mod} with formula \code{y ~ x1 + x2 + x3}, the predictor effect \code{p1 <- predictorEffects(mod, ~ x1)} is essentially equilavent to \code{p2 <- Effect("x1", mod)}.  When plotted, these objects may produce different graphs because \code{plot(p1)} will always put \code{x1} on the horizontal axis while \code{plot(p2)} uses a rule to determine the horizontal axis based on the characteristics of all the predictors, e.g., preferring continuous predictors over factors.

If \code{mod} has the formula \code{y ~ x1 + x2 + x3 + x1:x2}, then \code{p1 <- predictorEffects(mod, ~ x1)} is essentially equivalent to \code{p2 <- Effect(c("x1", "x2"), mod)}.  As in the last example, the plotted versions of these objects may differ because of different rules used to determine the predictor on the horizontal axis.

If \code{mod} has the formula \code{y ~ x1 + x2 + x3 + x1:x2 + x1:x3}, then \code{p1 <- predictorEffects(mod, ~ x1)} is essentially equilavent to \code{p2 <- Effect(c("x1", "x2", "x3"), mod)}.  Again, the plotted versions of these objects may differ because of the rules used to determine the horizontal axis.
}
\value{
\code{predictorEffect} returns an object of class \code{c("predictoreff", "eff")}.  The components of the object are described in the help for \code{\link{Effect}}; \code{predictorEffects} returns an object of class \code{"predictorefflist"}, which is a list whose elements are of class \code{c("predictoreff", "eff")}.
}
\references{
See \code{\link{Effect}}.
}
\author{
S. Weisberg, \email{sandy@umn.edu}
}



\seealso{
\code{\link{Effect}}, \code{\link{plot.predictoreff}}
}
\examples{
mod <- lm(prestige ~ type*(education + income) + women, Prestige)
plot(predictorEffect("income", mod))
plot(predictorEffects(mod, ~ education + income + women))

# svyglm() example (adapting an example from the survey package)
\donttest{
if (require(survey)){
  data(api)
  dstrat<-svydesign(id=~1, strata=~stype, weights=~pw,
    data=apistrat, fpc=~fpc)
  mod <- svyglm(sch.wide ~ ell + meals + mobility, design=dstrat,
    family=quasibinomial())
  plot(predictorEffects(mod),
    axes=list(y=list(lim=log(c(0.4, 0.99)/c(0.6, 0.01)),
      ticks=list(at=c(0.4, 0.75, 0.9, 0.95, 0.99)))))
}
}
}

\keyword{hplot}
\keyword{models}
