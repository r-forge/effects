\name{effect}
\alias{effect}
\alias{effect.default}
\alias{Effect}
\alias{Effect.default}
\alias{Effect.lm}
\alias{Effect.multinom}
\alias{Effect.merMod}
\alias{Effect.mlm}
\alias{Effect.poLCA}
\alias{Effect.polr}
\alias{Effect.svyglm}
\alias{allEffects}
\alias{allEffects.default}

\title{Functions For Constructing Effect Displays} 
\description{
  \code{Effect} and \code{effect} construct an \code{"eff"} object for a term (usually a high-order term) in a regression that models a response as a linear function of main effects and interactions of factors and covariates. These models include, among others, linear models (fit by \code{\link[stats]{lm}} and \code{\link[nlme]{gls}}), and generalized linear models (fit by \code{\link[stats]{glm}}), for which an \code{"eff"} object is created, and multinomial and proportional-odds logit models (fit respectively by \code{\link[nnet]{multinom}} and \code{\link[MASS]{polr}}), for which an \code{"effpoly"} object is created. The computed effect absorbs the lower-order terms marginal to the term in question, and averages over other terms in the model.  For multivariate linear models (of class \code{"mlm"}, fit by \code{\link[stats]{lm}}), the functions construct a list of \code{"eff"} objects, separately for the various response variables in the model. 

\code{effect} builds the required object by specifying explicitly a focal term like \code{"a:b"} for an \code{a} by \code{b} interaction.  \code{Effect} in contrast specifies the predictors in  a term, for example \code{c("a", "b")}, rather than the term itself.  \code{Effect} is consequently more flexible and robust than \code{effect}, and will succeed with some models for which \code{effect} fails.  The \code{effect} function works by constructing a call to \code{Effect} and continues to be included in \pkg{effects} so older code that uses it will not break.

The \code{Effect} and \code{effect} functions can also be used with many other models; see \code{\link{Effect.default}} and the \href{../doc/methods-supported-by-effects.pdf}{Regression Models Supported by the effects Package} vignette.

\code{allEffects} identifies all of the high-order terms in a model and returns a list of \code{"eff"} or \code{"effpoly"} objects (i.e., an object of class \code{"efflist"}).

For information on computing and displaying \emph{predictor effects}, see \code{\link{predictorEffect}} and \code{\link{plot.predictoreff}}.

For further information about plotting effects, see \code{\link{plot.eff}}.
}

\usage{
effect(term, mod, vcov.=vcov, ...)

\method{effect}{default}(term, mod, vcov.=vcov, ...)

Effect(focal.predictors, mod, ...)

\method{Effect}{lm}(focal.predictors, mod, xlevels=list(),
    fixed.predictors, vcov. = vcov, se=TRUE,
    residuals=FALSE, quantiles=seq(0.2, 0.8, by=0.2),
    x.var=NULL, transformation, ...,
    #legacy arguments:
    given.values, typical, offset, confint, confidence.level,
    partial.residuals)

\method{Effect}{multinom}(focal.predictors, mod,
    xlevels=list(), fixed.predictors,
    vcov. = vcov, se=TRUE, ...,
    #legacy arguments:
    confint, confidence.level, given.values, typical)

\method{Effect}{polr}(focal.predictors, mod,
    xlevels=list(), fixed.predictors,
    vcov.=vcov, se=TRUE, latent=FALSE, ...,
    #legacy arguments:
    confint, confidence.level, given.values, typical)

\method{Effect}{svyglm}(focal.predictors, mod, fixed.predictors, ...)

\method{Effect}{merMod}(focal.predictors, mod, ..., KR=FALSE)
    
\method{Effect}{poLCA}(focal.predictors, mod, ...)

\method{Effect}{mlm}(focal.predictors, mod, response, ...)

allEffects(mod, ...)

\method{allEffects}{default}(mod, ...)
}

\arguments{
  
  \item{term}{the quoted name of a term, usually, but not necessarily, a high-order term in the model. The term must be given exactly as it appears in the printed model, although either colons (\code{:}) or asterisks (\code{*}) may be used for interactions.  If \code{term} is NULL, the function returns the formula for the linear predictor.}
  
  \item{focal.predictors}{a character vector of one or more predictors in the model in any order.}
  
  \item{mod}{a regression model object.  If no specific method exists for the class of \code{mod}, \code{Effect.default} will be called.}
  
  \item{xlevels}{this argument is used to set the number of levels for any focal numeric predictor (that is predictors that are not factors, character variables, or logical variables, all of which are treated as factors).  If \code{xlevels=NULL}, then each numeric predictor is represented by five values over its range, equally spaced and then rounded to 'nice' numbers.  If \code{xlevels=n} is an integer, then each numeric predictor is represented by \code{n} equally spaced values rounded to 'nice' numbers.  
    
    More generally, \code{xlevels} can be a named list of values at which to set each numeric predictor.  For example, \code{xlevels=list(x1=c(2, 4.5, 7), x2=4)} would use the values 2, 4.5, and 7 for \code{x1}, use 4 equally spaced values for \code{x2}, and use the default for any other numeric predictors. 
    
    If partial residuals are computed, then the focal predictor that is to appear on the horizontal axis of an effect plot is evaluated at 100 equally spaced values along its full range, and, by default, other numeric predictors are evaluated at the quantiles specified in the \code{quantiles} argument, unless their values are given explicitly in \code{xlevels}.}
  
  \item{fixed.predictors}{an optional list of specifications affecting the values at which fixed predictors for an effect are set, potentially including:

    \describe{

      \item{given.values}{\code{given.values="default"} (which is, naturally, the default) specifies averaging over levels of a non-focal factor, weighting levels of the factor in proportion to sample size.
      
      \code{given.values="equal"} computes unweighted averages over the levels of non-focal factors.  
      
      For finer control, the user can also provide a named numeric vector of weights for particular columns of the model matrix that correspond to the regressors for the factor. 
      
Character and logical predictors are treated as factors. 

For example, for a factor \code{X} with three levels \code{a}, \code{b} and \code{c}, the regressors generated using the default \code{\link[stats]{contr.treatment}} parameterization for a factor will be named \code{Xb} and \code{Xc}, as the regressor for level \code{a} is excluded as the baseline level.  The specification \code{given.values=c(Xb=1/2, Xc=1/4)} would average over the levels of \code{X} with weight 1/2 for level \code{b}, 1/4 for \code{c}, and weight 1 = 1/2 - 1/4 = 1/4 for the baseline level \code{a}.  Setting \code{given.values=c(Xb=1)} would fix \code{X} at level \code{b}.
}

  	  \item{typical}{a function to be applied to the columns of the model matrix over which the effect is "averaged"; with the exception of the \code{"svyglm"} method, the default is \code{\link{mean}}. For\code{"svyglm"} objects, the default is to use the survey-design weighted mean.}

      \item{apply.typical.to.factors}{It generally doesn't make sense to apply typical values that aren't means (e.g., medians) to the columns of the model-matrix representing contrasts for factors. This value generally defaults to \code{FALSE} except for \code{"svyglm"} objects, for which the default is \code{TRUE}, using the the survey-design weighted mean.}

      \item{offset}{a function to be applied to the offset values (if there is an offset) in a linear or generalized linear model, or a mixed-effects model fit by \code{\link[lme4]{lmer}} or \code{\link[lme4]{glmer}}; or a numeric value, to which the offset will be set. The default is the \code{\link{mean}} function, and thus the offset will be set to its mean; in the case of \code{"svyglm"} objects, the default is to use the survey-design weighted mean. \emph{Note:} Only offsets defined by the \code{offset} argument to \code{\link[stats]{lm}}, \code{\link[stats]{glm}}, \code{\link[survey]{svyglm}}, \code{\link[lme4]{lmer}}, or \code{\link[lme4]{glmer}} will be handled correctly; use of the \code{offset} function in the model formula is not supported.}
  }
}

  \item{vcov.}{Effect methods generally use the matrix returned by \code{vcov(mod)} to compute standard errors and confidence bounds.  Alternatively, the user may specify the name of a function that returns a matrix of the same dimension and structure as the matrix returned by \code{vcov(mod)}.  For example, \code{vcov. = hccm} uses the \code{\link[car]{hccm}} function from the \pkg{car} package to use a heteroscedasticity corrected covariance matrix for a linear model in place of the standard covariance estimate. This argument can be set to equal matrix of the same size and structure as the matrix returned by \code{vcov(mod)}. For example, using \code{vcov. = vcov(Boot(mod))} uses \code{\link[car]{Boot}} from the \pkg{car} package to get a bootstrap estimate of the covariance matrix for linear, generalized linear, and possibly other modeling frameworks.}
  
  \item{se}{\code{TRUE} (the default), \code{FALSE}, or a list with any or all of the following elements, controlling whether and how standard errors and confidence limits are computed for the effects: 
  \describe{
  
  \item{compute}{(default \code{TRUE}) whether or not to compute standard errors and confidence limits.} 
  
  \item{level}{(default \code{0.95}) confidence level for confidence limits.}
  
  \item{type}{one of \code{"pointwise"} (the default), \code{"Scheffe"}, or \code{"scheffe"}, whether to compute confidence limits with specified coverage at each point for an effect or to compute limits for a Scheffe-type confidence envelope. For \code{mer}, \code{merMod}, and \code{lme} objects, the normal distribution is used to get confidence limits.}
  }
}

  \item{residuals}{if \code{TRUE}, residuals for a linear or generalized linear model will be computed and saved; if \code{FALSE} (the default), residuals are suppressed. If residuals are saved, partial residuals are computed when the effect is plotted: see \code{\link{plot.eff}} and the vignette \href{../doc/partial-residuals.pdf}{Effect Displays with Partial Residuals}. This argument may also be used for mixed-effects and some other models.}
  
  \item{quantiles}{quantiles at which to evaluate numeric focal predictors \emph{not} on the horizontal axis, used only when partial residuals are displayed; superseded if the \code{xlevels} argument gives specific values for a predictor.}
    
  \item{x.var}{the (quoted) name or index of the numeric predictor to define the horizontal axis of an effect plot for a linear or generalized linear model; the default is \code{NULL}, in which case the first numeric predictor in the effect will be used \emph{if} partial residuals are to be computed. This argument is intended to be used when \code{residuals} is \code{TRUE}; otherwise, the variable on the horizontal axis can be chosen when the effect object is plotted: see \code{\link{plot.eff}}.}
  
  \item{transformation}{for the \code{Effect.lm} method, an optional two-element list with \code{link} and \code{inverse} elements to transform the response (see examples); an alternative to use for graphs is to set the argument \code{axes = list(y = list(transformation = list(link = link-function, inverse = mean-function)))} (see \code{\link{plot.eff})}; the argument must be used for transforming the response in printed or summary output.}
  
  \item{latent}{if \code{TRUE}, effects in a proportional-odds logit model are computed on the scale of the latent response; if \code{FALSE}	(the default) effects are computed as individual-level probabilities and logits.}
  
  \item{x}{an object of class \code{"eff"}, \code{"effpoly"}, or \code{"efflatent"}.}
  
  \item{KR}{if \code{TRUE} and the \pkg{pbkrtest} package is installed, use the Kenward-Roger coefficient covariance matrix to compute effect standard errors for linear mixed models fit with \code{\link[lme4]{lmer}}; the default is \code{FALSE} because the computation can be time-consuming.}
  
  \item{response}{for an \code{"mlm"} object, a vector containing the (quoted) name(s) or indices of one or more response variable(s). The default is to use all responses in the model.}
    
  \item{...}{arguments to be passed down.}
  
  \item{confint, confidence.level, given.values, typical, offset, partial.residuals}{legacy arguments retained for backwards compatibility; if present, these arguments take precedence over the \code{level} element of the \code{confint} list argument and the \code{given.values}, \code{typical}, and \code{offset} elements of the \code{fixed.predictors} list argument; \code{confint} may be used in place of the \code{se} argument; \code{partial.residuals} may be used in place of the \code{residuals} argument. See \code{\link{LegacyArguments}} for details.}
}

\details{
  Normally, the functions to be used directly are \code{allEffects}, to return a list of high-order effects, and the generic \code{plot} function to plot the effects (see \code{\link{plot.efflist}}, \code{\link{plot.eff}}, and \code{\link{plot.effpoly}}). Alternatively, \code{Effect} can be used to vary a subset of predictors over their ranges, while other predictors are held to typical values.
  
Plotting methods for effect objects call the \code{\link[lattice]{xyplot}} (or in some cases, the \code{\link[lattice]{densityplot}}) function in the \pkg{lattice} package. Effects may also be printed (implicitly or explicitly via \code{print}) or summarized (using \code{summary}) (see \code{\link{print.efflist}}, \code{\link{summary.efflist}}, \code{\link{print.eff}}, \code{\link{summary.eff}}, \code{\link{print.effpoly}}, and \code{\link{summary.effpoly}}).

  If asked, the \code{effect} function will compute effects for terms that have higher-order relatives in the model, averaging over those terms (which rarely makes sense), or for terms that do not appear in the model but are higher-order relatives of terms that do. For example, for the model \code{Y ~ A*B + A*C + B*C}, one could compute the effect corresponding to the absent term \code{A:B:C}, which absorbs the constant, the \code{A}, \code{B}, and \code{C} main effects, and the three two-way interactions. In either of these cases, a warning is printed.

  See \code{\link{predictorEffects}} for an alternative paradigm for defining effects.
}

\value{

  For \code{"lm"}, \code{"glm"}, \code{"svyglm"}, \code{"lmerMod"}, \code{"glmerMod"}, and \code{"lme"}, model objects, \code{effect} and \code{Effect} return  an \code{"eff"} object, and for \code{"multinom"}, \code{"polr"}, \code{"clm"}, \code{"clmm"}, and \code{"clm2"} models, an \code{"effpoly"} object, with the components listed below. For an \code{"mlm"} object with one response specified, an \code{"eff"} object is returned, otherwise an \code{"efflist"} object is returned, containing one \code{"eff"} object for each \code{response}.
  
  \item{term}{the term to which the effect pertains.}
  
  \item{formula}{the complete model formula.}
  
  \item{response}{a character string giving the name of the response variable.}
  
  \item{y.levels}{(for \code{"effpoly"} objects) levels of the polytomous response variable.}
  
  \item{variables}{a list with information about each predictor, including its name, whether it is a factor, and its levels or values.}
  
  \item{fit}{(for \code{"eff"} objects) a one-column matrix of fitted values, representing the effect on the scale of the linear predictor; this is a raveled table, representing all combinations of predictor values.}
  
  \item{prob}{(for \code{"effpoly"} objects) a matrix giving fitted probabilities for the effect for the various levels of the the response (columns) and combinations of the focal predictors (rows).}
  
  \item{logit}{(for \code{"effpoly"} objects) a matrix giving fitted logits for the effect for the various levels of the the response (columns) and combinations of the focal predictors (rows).}
  
  \item{x}{a data frame, the columns of which are the predictors in the effect, and the rows of which give all combinations of values of these predictors.}
  
  \item{model.matrix}{the model matrix from which the effect was calculated.}
  
  \item{data}{a data frame with the data on which the fitted model was based.}
  
  \item{discrepancy}{the percentage discrepancy for the `safe' predictions of the original fit; should be very close to 0.  Note: except for \code{gls} models, this is now necessarily 0.}
  
  \item{offset}{value to which the offset is fixed; \code{0} if there is no offset.}
  
  \item{model}{(for \code{"effpoly"} objects) \code{"multinom"} or \code{"polr"}, as appropriate.}
  
  \item{vcov}{(for \code{"eff"} objects) a covariance matrix for the effect, on the scale of the linear predictor.}
  
  \item{se}{(for \code{"eff"} objects) a vector of standard errors for the effect, on the scale of the linear predictor.}
  
  \item{se.prob, se.logit}{(for \code{"effpoly"} objects) matrices of standard errors for the effect, on the probability and logit scales.}
  
  \item{lower, upper}{(for \code{"eff"} objects) one-column matrices of confidence limits, on the scale of the linear predictor.}
  
  \item{lower.prob, upper.prob, lower.logit, upper.logit}{(for \code{"effpoly"} objects) matrices of confidence limits for the fitted logits and probabilities; the latter are computed by transforming the former.}
  
  \item{confidence.level}{for the confidence limits.}
  
  \item{transformation}{(for \code{"eff"} objects) a two-element list, with element \code{link} giving the link function, and element \code{inverse} giving the inverse-link (mean) function; may be set directly via the \code{transformation} argument to \code{Effect.lm} or inferred from the model.}
  
  \item{residuals}{(working) residuals for linear or generalized linear models (and some similar models), to be used by \code{\link{plot.eff}} to compute and plot partial residuals.}
  
  \item{x.var}{the name of the predictor to appear on the horizontal axis of an effect plot made from the  returned object; will usually be \code{NULL} if partial residuals aren't computed.}
  
  \item{family}{for a \code{"glm"} model, the name of the distributional family of the model; for an \code{"lm"} model, this is \code{"gaussian"}; otherwise \code{NULL}. The \code{family} controls how partial residuals are smoothed in plots.}
  
  \item{link}{the value returned by \code{family(mod)}.  Down-stream methods may need the link, inverse link and derivative functions.}

  \code{allEffects} returns an \code{"efflist"} object, a list of \code{"eff"} or \code{"effpoly"} objects corresponding to the high-order terms of the model.

If \code{mod} is of class \code{"poLCA"} (from the \pkg{poLCA} package), representing a polytomous latent class model, effects are computed for the predictors given the estimated latent classes. The result is of class \code{"eff"} if the latent class model has 2 categories and of class \code{"effpoly"} with more than 2 categories.
}

\section{Warnings and Limitations}{

  The \code{Effect} function handles factors and covariates differently, and is likely to be confused if one is changed to the other in a model formula. Consequently, formulas that include calls to \code{as.factor}, \code{factor}, or \code{numeric} (as, e.g., in \code{y ~ as.factor(income)}) will cause errors. Instead, create the modified variables outside of the model formula (e.g., \code{fincome <- as.factor(income)}) and use these in the model formula.

The \code{effect} function doesn't work with factors that have colons in level names (e.g., \code{"level:A"}); the \code{effect} function will confuse the colons with interactions; rename levels to remove or replace the colons (e.g., \code{"level.A"}). Level names with colons are perfectly fine for use with \code{Effect}.

The functions in the \pkg{effects} package work properly with predictors that are numeric variables, factors, character variables, or logical variables; consequently, e.g., convert dates to numeric. Character predictors and logical predictors are treated as factors, the latter with "levels" \code{"FALSE"} and \code{"TRUE"}.

  Empty cells in crossed-factors are now permitted for \code{"lm"}, \code{"glm"}, and \code{"multinom"} models. For \code{"multinom"} models with two or more crossed factors with an empty cell, stacked area plots apparently do not work because of a bug in the \code{\link[lattice]{barchart}} function in the \pkg{lattice} package.  However, the default line plots do work.

  Offsets in linear and generalized linear models are supported, as are offsets in mixed models fit by \code{lmer} or \code{glmer}, but must be supplied through the \code{offset} argument to \code{lm}, \code{glm}, \code{lmer} or \code{glmer}; offsets supplied via calls to the \code{offset} function on the right-hand side of the model formula are not supported.

  Fitting ordinal mixed models using \code{\link[ordinal]{clmm}} or \code{\link[ordinal:clmmOld]{clmm2}} permits many options, including a variety of link functions, scale functions, nominal regressors, and various methods for setting thresholds.  Effects are currently generated only for the default values of the arguments \code{scale}, \code{nominal}, \code{link}, and \code{threshold}, which is equivalent to fitting an ordinal-response mixed-effects model with a logit link. \code{Effect} can also be used with objects created by \code{\link[ordinal]{clm}} or \code{\link[ordinal:clmOld]{clm2}}, fitting ordinal response models with the same links permitted by \code{\link[MASS]{polr}} in the \pkg{MASS} package, with no random effects, and with results similar to those from \code{\link[MASS]{polr}}.

  Calling any of these functions from within a user-written function may result in errors due to R's scoping rules.  See the vignette \code{embedding.pdf} in the \pkg{car} package for a solution to this problem.
  }

\references{
    Fox, J. (1987).
    Effect displays for generalized linear models.
    \emph{Sociological Methodology}
    \bold{17}, 347--361.

    Fox, J. (2003)
    Effect displays in R for generalised linear models.
    \emph{Journal of Statistical Software}
    \bold{8:15}, 1--27, \doi{10.18637/jss.v008.i15}.

  	Fox, J. and R. Andersen (2006).
   	Effect displays for multinomial and proportional-odds logit models.
  	\emph{Sociological Methodology}
  	\bold{36}, 225--255.

  	Fox, J. and J. Hong (2009).
    Effect displays in R for multinomial and proportional-odds logit models:?
    	    Extensions to the effects package.
    \emph{Journal of Statistical Software}
    \bold{32:1}, 1--24, \doi{10.18637/jss.v032.i01}.
    
    Fox, J. and S. Weisberg (2019).  \emph{An R Companion to Applied Regression, third edition},
    Thousand Oaks: Sage.

    Fox, J. and S. Weisberg (2018).
    Visualizing Fit and Lack of Fit in Complex Regression Models
        with Predictor Effect Plots with Partial Residuals.
    \emph{Journal of Statistical Software}
    \bold{87:9}, 1--27, \doi{10.18637/jss.v087.i09}.

    Hastie, T. J. (1992).
    Generalized additive models.
    In Chambers, J. M., and Hastie, T. J. (eds.)
    \emph{Statistical Models in S}, Wadsworth.

    Weisberg, S. (2014).
    \emph{Applied Linear Regression}, 4th edition, Wiley,
    \url{http://z.umn.edu/alr4ed}.
    }



\author{John Fox \email{jfox@mcmaster.ca}, Sanford Weisberg \email{sandy@umn.edu}
and Jangman Hong.}


\seealso{\code{\link{LegacyArguments}}. For information on printing, summarizing, and plotting effects:
    \code{\link{print.eff}}, \code{\link{summary.eff}}, \code{\link{plot.eff}},
    \code{\link{print.summary.eff}},
    \code{\link{print.effpoly}}, \code{\link{summary.effpoly}}, \code{\link{plot.effpoly}},
    \code{\link{print.efflist}}, \code{\link{summary.efflist}},
    \code{\link{plot.efflist}}, \code{\link[lattice]{xyplot}},
    \code{\link[lattice]{densityplot}}, and the \href{../doc/partial-residuals.pdf}{Effect Displays with Partial Residuals} and \href{../doc/methods-supported-by-effects.pdf}{Regression Models Supported by the effects Package} vignettes.}

\examples{
mod.cowles <- glm(volunteer ~ sex + neuroticism*extraversion,
                  data=Cowles, family=binomial)
eff.cowles <- allEffects(mod.cowles, xlevels=list(extraversion=seq(0, 24, 6)),
                   fixed.predictors=list(given.values=c(sexmale=0.5)))
eff.cowles
as.data.frame(eff.cowles[[2]])

\donttest{
# the following are equivalent:
eff.ne <- effect("neuroticism*extraversion", mod.cowles)
Eff.ne <- Effect(c("neuroticism", "extraversion"), mod.cowles)
all.equal(eff.ne$fit, Eff.ne$fit)

plot(eff.cowles, 'sex', axes=list(y=list(lab="Prob(Volunteer)")))

plot(eff.cowles, 'neuroticism:extraversion',
     axes=list(y=list(lab="Prob(Volunteer)",
        ticks=list(at=c(.1,.25,.5,.75,.9)))))

plot(Effect(c("neuroticism", "extraversion"), mod.cowles,
            se=list(type="Scheffe"),
            xlevels=list(extraversion=seq(0, 24, 6)),
            fixed.predictors=list(given.values=c(sexmale=0.5))),
     axes=list(y=list(lab="Prob(Volunteer)",
        ticks=list(at=c(.1,.25,.5,.75,.9)))))

plot(eff.cowles, 'neuroticism:extraversion', lines=list(multiline=TRUE),
     axes=list(y=list(lab="Prob(Volunteer)")))

plot(effect('sex:neuroticism:extraversion', mod.cowles,
            xlevels=list(extraversion=seq(0, 24, 6))),
     lines=list(multiline=TRUE))
}

# a nested model:

mod <- lm(log(prestige) ~ income:type + education, data=Prestige)

plot(Effect(c("income", "type"), mod, transformation=list(link=log, inverse=exp)),
     axes=list(y=list(lab="prestige")))


if (require(nnet)){
    mod.beps <- multinom(vote ~ age + gender + economic.cond.national +
                             economic.cond.household + Blair + Hague + Kennedy +
                             Europe*political.knowledge, data=BEPS)
    \donttest{
    plot(effect("Europe*political.knowledge", mod.beps,
                xlevels=list(political.knowledge=0:3)))
    }

    plot(Effect(c("Europe", "political.knowledge"), mod.beps,
                xlevels=list(Europe=1:11, political.knowledge=0:3),
                fixed.predictors=list(given.values=c(gendermale=0.5))),
         lines=list(col=c("blue", "red", "orange")),
         axes=list(x=list(rug=FALSE), y=list(style="stacked")))

    \donttest{
    plot(effect("Europe*political.knowledge", mod.beps, # equivalent
                xlevels=list(Europe=1:11, political.knowledge=0:3),
                fixed.predictors=list(given.values=c(gendermale=0.5))),
         lines=list(col=c("blue", "red", "orange")),
         axes=list(x=list(rug=FALSE), y=list(style="stacked")))
    }
}

if (require(MASS)){
    mod.wvs <- polr(poverty ~ gender + religion + degree + country*poly(age,3),
                    data=WVS)
    \donttest{
    plot(effect("country*poly(age, 3)", mod.wvs))
    }

    plot(Effect(c("country", "age"), mod.wvs),
         axes=list(y=list(style="stacked")))

    \donttest{
    plot(effect("country*poly(age, 3)", mod.wvs),
         axes=list(y=list(style="stacked"))) # equivalent

    plot(effect("country*poly(age, 3)", latent=TRUE, mod.wvs))
    plot(effect("country*poly(age, 3)", latent=TRUE, mod.wvs,
         se=list(type="scheffe"))) # Scheffe-type confidence envelopes
    }
}


mod.pres <- lm(prestige ~ log(income, 10) + poly(education, 3) + poly(women, 2),
               data=Prestige)
eff.pres <- allEffects(mod.pres, xlevels=50)
plot(eff.pres)
plot(eff.pres[1],
     axes=list(x=list(income=list(
             transform=list(trans=log10, inverse=function(x) 10^x),
             ticks=list(at=c(1000, 2000, 5000, 10000, 20000))
    ))))
\donttest{
# linear model with log-response and log-predictor
# to illustrate transforming axes and setting tick labels
mod.pres1 <- lm(log(prestige) ~ log(income) + poly(education, 3) + poly(women, 2),
                data=Prestige)
# effect of the log-predictor
eff.log <- Effect("income", mod.pres1)
# effect of the log-predictor transformed to the arithmetic scale
eff.trans <- Effect("income", mod.pres1, transformation=list(link=log, inverse=exp))
#variations:
# y-axis:  scale is log, tick labels are log
# x-axis:  scale is arithmetic, tick labels are arithmetic
plot(eff.log)

# y-axis:  scale is log, tick labels are log
# x-axis:  scale is log, tick labels are arithmetic
plot(eff.log, axes=list(x=list(income=list(
    transform=list(trans=log, inverse=exp),
    ticks=list(at=c(5000, 10000, 20000)),
    lab="income, log-scale"))))

# y-axis:  scale is log, tick labels are arithmetic
# x-axis:  scale is arithmetic, tick labels are arithmetic
plot(eff.trans, axes=list(y=list(lab="prestige")))

# y-axis:  scale is arithmetic, tick labels are arithmetic
# x-axis:  scale is arithmetic, tick labels are arithmetic
plot(eff.trans, axes=list(y=list(type="response", lab="prestige")))

# y-axis:  scale is log, tick labels are arithmetic
# x-axis:  scale is log, tick labels are arithmetic
plot(eff.trans, axes=list(
       x=list(income=list(
            transform=list(trans=log, inverse=exp),
            ticks=list(at=c(1000, 2000, 5000, 10000, 20000)),
            lab="income, log-scale")),
        y=list(lab="prestige, log-scale")),
     main="Both response and X in log-scale")

# y-axis:  scale is arithmetic, tick labels are arithmetic
# x-axis:  scale is log, tick labels are arithmetic
plot(eff.trans, axes=list(
        x=list(
            income=list(transform=list(trans=log, inverse=exp),
                        ticks=list(at=c(1000, 2000, 5000, 10000, 20000)),
                        lab="income, log-scale")),
        y=list(type="response", lab="prestige")))
}

if (require(nlme)){ # for gls()
    mod.hart <- gls(fconvict ~ mconvict + tfr + partic + degrees, data=Hartnagel,
                    correlation=corARMA(p=2, q=0), method="ML")
    plot(allEffects(mod.hart))
    detach(package:nlme)
}

if (require(lme4)){
    data(cake, package="lme4")
    fm1 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake,
                REML = FALSE)
    plot(Effect(c("recipe", "temperature"), fm1))
    \donttest{
    plot(effect("recipe:temperature", fm1),
         axes=list(grid=TRUE)) # equivalent (plus grid)
    }
    if (any(grepl("pbkrtest", search()))) detach(package:pbkrtest)
    detach(package:lme4)
}

\donttest{
if (require(nlme) && length(find.package("lme4", quiet=TRUE)) > 0){
    data(cake, package="lme4")
    cake$rep <- with(cake, paste( as.character(recipe), as.character(replicate), sep=""))
    fm2 <- lme(angle ~ recipe * temperature, data=cake,
               random = ~ 1 | rep, method="ML")
    plot(Effect(c("recipe", "temperature"), fm2))
    plot(effect("recipe:temperature", fm2),
         axes=list(grid=TRUE))  # equivalent (plus grid)
    }
    detach(package:nlme)
}

\donttest{
if (require(poLCA)){
    data(election)
    f2a <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
                 MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~PARTY*AGE
    nes2a <- poLCA(f2a,election,nclass=3,nrep=5)
    plot(Effect(c("PARTY", "AGE"), nes2a),
         axes=list(y=list(style="stacked")))
}
}

# mlm example
if (require(heplots)) {
    data(NLSY, package="heplots")
    mod <- lm(cbind(read,math) ~ income+educ, data=NLSY)
    eff.inc <- Effect("income", mod)
    plot(eff.inc)
    eff.edu <- Effect("educ", mod)
    plot(eff.edu, axes=list(x=list(rug=FALSE), grid=TRUE))
    \donttest{
    plot(Effect("educ", mod, response="read"))
    }
    detach(package:heplots)
}

# svyglm() example (adapting an example from the survey package)

\donttest{
if (require(survey)){
  data("api")
  dstrat<-svydesign(id=~1, strata=~stype, weights=~pw,
    data=apistrat, fpc=~fpc)
  mod <- svyglm(sch.wide ~ ell + meals + mobility, design=dstrat,
    family=quasibinomial())
  plot(allEffects(mod),
    axes=list(y=list(lim=log(c(0.4, 0.99)/c(0.6, 0.01)),
      ticks=list(at=c(0.4, 0.75, 0.9, 0.95, 0.99)))))
}
}

# component + residual plot examples

\donttest{
Prestige$type <- factor(Prestige$type, levels=c("bc", "wc", "prof"))

mod.prestige.1 <- lm(prestige ~ income + education, data=Prestige)
plot(allEffects(mod.prestige.1, residuals=TRUE)) # standard C+R plots
plot(allEffects(mod.prestige.1, residuals=TRUE,
    se=list(type="scheffe"))) # with Scheffe-type confidence bands

mod.prestige.2 <- lm(prestige ~ type*(income + education), data=Prestige)
plot(allEffects(mod.prestige.2, residuals=TRUE))

mod.prestige.3 <- lm(prestige ~ type + income*education, data=Prestige)
plot(Effect(c("income", "education"), mod.prestige.3, residuals=TRUE),
     partial.residuals=list(span=1))
}

#  artificial data

set.seed(12345)
x1 <- runif(500, -75, 100)
x2 <- runif(500, -75, 100)
y <- 10 + 5*x1 + 5*x2 + x1^2 + x2^2 + x1*x2 + rnorm(500, 0, 1e3)
Data <- data.frame(y, x1, x2)
mod.1 <- lm(y ~ poly(x1, x2, degree=2, raw=TRUE), data=Data)
# raw=TRUE necessary for safe prediction
mod.2 <- lm(y ~ x1*x2, data=Data)
mod.3 <- lm(y ~ x1 + x2, data=Data)

plot(Effect(c("x1", "x2"), mod.1, residuals=TRUE)) # correct model
plot(Effect(c("x1", "x2"), mod.2, residuals=TRUE)) # wrong model
plot(Effect(c("x1", "x2"), mod.3, residuals=TRUE)) # wrong model
}

\keyword{hplot}
\keyword{models}
