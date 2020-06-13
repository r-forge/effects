# 12/11/2017:  S. Weisberg.  This file contains all the Effect methods that call
# Effect.default.  Excluded are Effect.lm, Effect.polr, and Effect.multinom, 
# and for now Effect.svyglm.
# 06/08/2018: rewrote method for betareg, removing the 'link' argument from sources
# 11/28/2018:  modified Effect.gls to ignore the weights argument by 
#              deleting it from sources$call.
# 11/30/2018: fixed bug in Effect.merMod() specifying fam$family explicitly.
# 7/5/2019:  clm clm2 and clmm were not passing the estimated threshholds to polr
# 3/22/2020:  added Effect.glmmPQL (from MASS package)
# 4/27/2020:  require 'insight' package for find_formula and get_coefficients
#             so formula and coefficients are generally not needed
# 2020-06-13: fix typo (omitted ') in an error message



# lme method
#Effect.lme <- function(focal.predictors, mod, ...){
#  args <- list()
#  Effect.default(focal.predictors, mod, ..., sources=args)
#}



# new gls method
Effect.gls <- function(focal.predictors, mod, ...){
  cl <- mod$call
  cl$weights <- NULL
  args <- list(
    call = cl)
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# new glmmPQL method 3/22/2020
Effect.glmmPQL <- function(focal.predictors, mod, ...){
  args <- list(
    family = mod$family)
  Effect.default(focal.predictors, mod, ..., sources=args)
}

Effect.merMod <- function(focal.predictors, mod, ..., KR=FALSE){
  if (KR && !requireNamespace("pbkrtest", quietly=TRUE)){
    KR <- FALSE
    warning("pbkrtest is not available, KR set to FALSE")}
  fam <- family(mod)
  args <- list(
    family=fam,
    vcov = if (fam$family == "gaussian" && fam$link == "identity" && KR)
      as.matrix(pbkrtest::vcovAdj(mod)) else insight::get_varcov(mod))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# rlmer in robustlmm package
Effect.rlmerMod <- function(focal.predictors, mod, ...){
  args <- list(
    family=family(mod))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# clm in the ordinal package. clm is not supported by insight package
Effect.clm <- function(focal.predictors, mod, ...){
  if (requireNamespace("MASS", quietly=TRUE)){
    polr <- MASS::polr} else stop("MASS package is required")
  polr.methods <- c("logistic", "probit", "loglog", 
                    "cloglog", "cauchit")
  method <- mod$link
  if(method == "logit") method <- "logistic"
  if(!(method %in% polr.methods)) 
    stop("'link' must be a 'method' supported by polr; see help(polr)")
  if(mod$threshold != "flexible") 
    stop("Effects only supports the 'flexible' threshold")
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  args <- list(
    type = "polr",
    coefficients = mod$beta,
    zeta = mod$alpha,  
    method=method,
    vcov = as.matrix(vcov(mod)[or, or]))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# clm2, this is supported by insight package
Effect.clm2 <- function(focal.predictors, mod, ...){
  if (requireNamespace("MASS", quietly=TRUE)){
      polr <- MASS::polr}
  polr.methods <- c("logistic", "probit", "loglog", 
                    "cloglog", "cauchit")
  method <- mod$link
  if(!(method %in% polr.methods)) 
    stop("'link' must be a 'method' supported by polr; see help(polr)")
  if(is.null(mod$Hessian)){
     message("\nRe-fitting to get Hessian\n")
     mod <- update(mod, Hess=TRUE)}
  if(mod$threshold != "flexible") 
    stop("Effects only supports the flexible threshold")
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  args <- list(
    type = "polr",
    formula = mod$call$location,
    coefficients = mod$beta,  
    zeta = mod$Theta,
    method=method,
    vcov = as.matrix(vcov(mod)[or, or]))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

#clmm in ordinal package
Effect.clmm <- function(focal.predictors, mod, ...){ 
  if (requireNamespace("MASS", quietly=TRUE)){
    polr <- MASS::polr}
  else stop("The MASS package must be installed")
  polr.methods <- c("logistic", "probit", "loglog", 
                    "cloglog", "cauchit")
  method <- mod$link
  if(method == "logit") method <- "logistic"
  if(!(method %in% polr.methods)) 
    stop("'link' must be a 'method' supported by polr; see help(polr)")
  if(is.null(mod$Hessian)){
    message("\nRe-fitting to get Hessian\n")
    mod <- update(mod, Hess=TRUE)}
  if(mod$threshold != "flexible") 
    stop("Only threshold='flexible' is supported by effects")
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  Vcov <- as.matrix(vcov(mod)[or, or])
  args <- list(
    type = "polr",
    formula = insight::find_formula(mod)$conditional,
    coefficients = mod$beta,
    zeta=mod$alpha,
    method=method,
    vcov = as.matrix(Vcov))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# betareg from the betareg package
Effect.betareg <- function(focal.predictors, mod, ...){
  coef <- mod$coefficients$mean
  vco <- vcov(mod)[1:length(coef), 1:length(coef)]
# betareg uses beta errors with mean link given in mod$link$mean.  
# Construct a family based on the binomial() family
  fam <- binomial(link=mod$link$mean)
# adjust the varince function to account for beta variance
  fam$variance <- function(mu){
    f0 <- function(mu, eta) (1-mu)*mu/(1+eta)
    do.call("f0", list(mu, mod$coefficient$precision))}
# adjust initialize
  fam$initialize <- expression({mustart <- y})
# collect arguments
  args <- list(
    call = mod$call,
    formula = formula(mod),
    family=fam,
    coefficients = coef,
    vcov = vco)
  Effect.default(focal.predictors, mod, ..., sources=args)
}

