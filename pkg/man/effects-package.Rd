\name{effects-package}
\Rdversion{1.1}
\alias{effects-package}
\alias{effects}
\docType{package}
\title{
Effect Displays for Linear, Generalized Linear, and Other Models
}
\description{
Graphical and tabular effect displays, e.g., of interactions, for 
various statistical models with linear predictors.
}
\details{
\tabular{ll}{
Package: \tab effects\cr
Version: \tab 4.2-2\cr
Date: \tab 2022-02-16\cr
Depends: \tab R (>= 3.5.0), carData\cr
Suggests: \tab pbkrtest (>= 0.4-4), nlme, MASS, poLCA, heplots, splines, ordinal, car, knitr, betareg, alr4, robustlmm\cr
Imports: \tab lme4, nnet, lattice, grid, colorspace, graphics, grDevices, stats, survey, utils, estimability, insight\cr
LazyLoad: \tab yes\cr
License: \tab GPL (>= 2)\cr
URL: \tab https://www.r-project.org, https://socialsciences.mcmaster.ca/jfox/\cr
}

This package creates effect displays for various kinds of models, as partly explained in the references.
Typical usage is \code{plot(allEffects(model))} or \code{plot(predictorEffects(model))},
where \code{model} is an appropriate fitted-model object.
Additional arguments to \code{\link{allEffects}}, \code{\link{predictorEffects}} and \code{\link[graphics]{plot}} can be used to customize the resulting
displays. The function \code{\link{effect}} can be employed to produce an effect display for a 
particular term in the model, or to which terms in the model are marginal. The function \code{\link{predictorEffect}} can be used
to construct an effect display for a particularly predictor. The function \code{\link{Effect}}
may similarly be used to produce an effect display for any combination of predictors. In any of the cases, use \code{plot}
to graph the resulting effect object. For linear and generalized linear models it is also possible to plot 
partial residuals to obtain (multidimensional) component+residual plots.
See \code{?effect}, \code{?Effect}, \code{?predictorEffect}, and \code{?plot.eff} for details.
}

\author{
John Fox, Sanford Weisberg, Brad Price, Michael Friendly, Jangman Hong, 
Robert Anderson, David Firth, Steve Taylor, and the R Core Team. 

Maintainer: John Fox <jfox@mcmaster.ca>
}

\references{
   Fox, J. and S. Weisberg (2019)
   \emph{An R Companion to Applied Regression, Third Edition}
   Sage Publications.
   
   Fox, J. (1987)
    Effect displays for generalized linear models.
    \emph{Sociological Methodology}
    \bold{17}, 347--361.
    
    Fox, J. (2003)
    Effect displays in R for generalised linear models.
    \emph{Journal of Statistical Software}
    \bold{8:15}, 1--27, \doi{10.18637/jss.v008.i15}.
    
    Fox, J. and R. Andersen (2006)
    Effect displays  for multinomial and proportional-odds logit models.
    \emph{Sociological Methodology} 
    \bold{36}, 225--255.
    
    Fox, J. and J. Hong (2009).
    Effect displays in R for multinomial and proportional-odds logit models: 
        Extensions to the effects package.
    \emph{Journal of Statistical Software}
    \bold{32:1}, 1--24, \doi{10.18637/jss.v032.i01}.
    
    Fox, J. and S. Weisberg (2018).
    Visualizing Fit and Lack of Fit in Complex Regression Models:
        Effect Plots with Partial Residuals.
    \emph{Journal of Statistical Software}
    \bold{87:9}, 1--27, \doi{10.18637/jss.v087.i09}.
}


\keyword{ package }

