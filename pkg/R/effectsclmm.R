# 2014-12-11 Effects plots for ordinal and ordinal mixed models from the 'ordinal' package
# 2014-12-11 effect.clm built from effect.mer as modified 2014-12-07,  by S. Weisberg
# 2015-06-10: requireNamespace("MASS") rather than require("MASS")
# 2016-02-12: added support for clmm and clm objects from 'ordinal' S. Weisberg
# 2016-08-16: added ... argument to effect() and Effect() methods. J. Fox
# 2017-03-28: fixed bug to allow data=m S. Weisberg
# 2017-11-17: allow links with clm rather than just logit.  S. Weisberg
# 2017-11-17: allow clm clm2 clmm to work with predictorEffects.  S. Weisberg
# 2017-11-17: add predictorEffect(s) methods for ordinal functions.  S. Weisberg

###
###  clm2
###
clm2.to.polr <- function(mod) {
  if (requireNamespace("MASS", quietly=TRUE)){
    polr <- MASS::polr
  }
  else stop("The MASS package is needed for this function")
  cl <- mod$call
  present <- match(c("scale", "nominal", "link", "threshold"), names(cl), 0L)
  if(any(present != 0)) {
    if(present[3] != 0){if(!(cl$link %in% c("logistic", "probit", "loglog", "cloglog", "cauchit")))
      stop("'link' must be logistic, probit, loglog, cloglog or cauchit to use with effects")}
#    if(present[3] != 0){if(cl$link != "logistic") stop("'link' must be 'logisitic' for use with effects")}
#    if(present[4] != 0){if(cl$threshold != "flexible") stop("'threshold' must be 'flexible' for use with effects")}
#    if(present[1] != 0){if(!is.null(cl$scale)) stop("'scale' must be NULL for use with effects")}
#    if(present[2] != 0){if(!is.null(cl$nominal)) stop("'nominal' must be NULL for use with effects")}
  }
  if(is.null(mod$Hessian)){
    message("\nRe-fitting to get Hessian\n")
    mod <- update(mod, Hess=TRUE)
  }
  cl$formula <- cl$location
  cl$method <- cl$link
  .m <- match(c("formula", "data", "subset","weights",
                "na.action",  "contrasts", "method"), names(cl), 0L)
  cl <- cl[c(1L, .m)]
  cl$start <- c(mod$beta, mod$Theta)
  cl[[1L]] <- as.name("polr")
  cl$control <- list(maxit=1)
  mod2 <- eval(cl)
  mod2$coefficients <- mod$beta
  # get vcov
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  mod2$vcov <- as.matrix(vcov(mod)[or, or])
  class(mod2) <- c("fakeclm2", class(mod2))
  mod2
}

#method for 'fakeglm' objects. Do not export
vcov.fakeclm2 <- function(object, ...) object$vcov

###
###   clmm
###
clmm.to.polr <- function(mod) {
  if (requireNamespace("MASS", quietly=TRUE)){
    polr <- MASS::polr
  }
  else stop("The MASS package is needed for this function")
  cl <- mod$call
#  present <- match(c("scale", "nominal", "link", "threshold"), names(cl), 0L)
#  if(any(present != 0)) {
#    if(present[3] != 0){if(cl$link != "logit") stop("'link' must be 'logit' for use with effects")}
#    if(present[4] != 0){if(cl$threshold != "flexible") stop("'threshold' must be 'flexible' for use with effects")}
#    if(present[1] != 0){if(!is.null(cl$scale)) stop("'scale' must be NULL for use with effects")}
#    if(present[2] != 0){if(!is.null(cl$nominal)) stop("'nominal' must be NULL for use with effects")}
#  }
  if(is.null(mod$Hessian)){
    message("\nRe-fitting to get Hessian\n")
    mod <- update(mod, Hess=TRUE)
  }
  cl$formula <- fixmod(as.formula(mod$formula))  # changed from clm2
  cl$method <- cl$link
  if(!is.null(cl$method)) {if(cl$method=="logit") cl$method="logistic"} 
  .m <- match(c("formula", "data", "subset","weights",
                "na.action",  "contrasts", "method"), names(cl), 0L)
  cl <- cl[c(1L, .m)]
  cl$start <- c(mod$beta, mod$Theta)
  cl$control <- list(maxit=1000)   ##########
  cl[[1L]] <- as.name("polr")
  cl$control <- list(maxit=1)
  mod2 <- eval(cl)
  mod2$coefficients <- mod$beta
  # get vcov
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  mod2$vcov <- as.matrix(vcov(mod)[or, or])
  class(mod2) <- c("fakeclmm", class(mod2))
  mod2
}

#method for 'fakeglm' objects. Do not export
vcov.fakeclmm <- function(object, ...) object$vcov


###
### clm
###
clm.to.polr <- function(mod) {
  if (requireNamespace("MASS", quietly=TRUE)){
    polr <- MASS::polr
  }
  else stop("The MASS package is needed for this function")
  cl <- mod$call
#  present <- match(c("scale", "nominal", "threshold"), names(cl), 0L)
#  if(any(present != 0)) {
#    if(present[3] != 0){if(cl$link != "logit") stop("'link' must be 'logit' for use with effects")}
#    if(present[4] != 0){if(cl$threshold != "flexible") stop("'threshold' must be 'flexible' for use with effects")}
#    if(present[1] != 0){if(!is.null(cl$scale)) stop("'scale' must be NULL for use with effects")}
#    if(present[2] != 0){if(!is.null(cl$nominal)) stop("'nominal' must be NULL for use with effects")}
#  }
  if(is.null(mod$Hessian)){
    message("\nRe-fitting to get Hessian\n")
    mod <- update(mod, Hess=TRUE)
  }
  cl$method <- cl$link
  .m <- match(c("formula", "data", "subset","weights",
                "na.action",  "contrasts", "method"), names(cl), 0L)
  cl <- cl[c(1L, .m)]
  cl$start <- c(mod$beta, mod$Theta)
  cl$control <- list(maxit=1)
  cl[[1L]] <- as.name("polr")
  mod2 <- eval(cl)
  mod2$coefficients <- mod$beta
  # get vcov
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  mod2$vcov <- as.matrix(vcov(mod)[or, or])
  class(mod2) <- c("fakeclm", class(mod2))
  mod2
}

#method for 'fakeglm' objects. Do not export
vcov.fakeclm <- function(object, ...) object$vcov

#The next three functions should be exported

effect.clm <- function(term, mod, ...) {
  effect(term, clm.to.polr(mod), ...)
}

allEffects.clm <- function(mod, ...){
  allEffects(clm.to.polr(mod), ...)
}

Effect.clm <- function(focal.predictors, mod, ...){
  Effect(focal.predictors, clm.to.polr(mod), ...)
}

#The next three functions should be exported

effect.clm2 <- function(term, mod, ...) {
  effect(term, clm2.to.polr(mod), ...)
}

allEffects.clm2 <- function(mod, ...){
  allEffects(clm.to.polr(mod), ...)
}

Effect.clm2 <- function(focal.predictors, mod, ...){
  Effect(focal.predictors, clm2.to.polr(mod), ...)
}

#The next three functions should be exported

effect.clmm <- function(term, mod, ...) {
  effect(term, clmm.to.polr(mod), ...)
}

allEffects.clmm <- function(mod, ...){
  allEffects(clmm.to.polr(mod), ...)
}

Effect.clmm <- function(focal.predictors, mod, ...){
  Effect(focal.predictors, clmm.to.polr(mod), ...)
}

predictorEffects.clm <- function(mod, predictors = ~., ...){
  predictorEffects(clm.to.polr(mod), predictors=predictors, ...)
}

predictorEffect.clm2 <- function(predictor, mod, xlevels=list(), ...){
  predictorEffect(predictor, clm2.to.polr(mod), xlevels=xlevels, ...)
}

predictorEffect.clmm <- function(predictor, mod, xlevels=list(), ...){
  predictorEffect(predictor, clmm.to.polr(mod), xlevels=xlevels, ...)
}

predictorEffects.clm2 <- function(mod, predictors = ~., ...){
  predictorEffects(clm2.to.polr(mod), predictors=predictors, ...)
}

predictorEffects.clmm <- function(mod, predictors = ~., ...){
  predictorEffects(clmm.to.polr(mod), predictors=predictors, ...)
}
