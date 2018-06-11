# 12/11/2017:  S. Weisberg.  This file contains all the Effect methods that call
# Effect.default.  Excluded are Effect.lm, Effect.polr, and Effect.multinom, 
# and for now Effect.svyglm.
# 06/08/2018: rewrote method for betareg, removing the 'link' argument from sources

# new lme method
Effect.lme <- function(focal.predictors, mod, ...){
  args <- list(
    call = mod$call,
    formula = mod$call$fixed,
    coefficients = mod$coefficients$fixed,
    vcov = mod$varFixed)
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# new gls method
Effect.gls <- function(focal.predictors, mod, ...){
  args <- list(
    call = mod$call,
    formula = formula(mod),
    coefficients = coef(mod),
    vcov = as.matrix(vcov(mod)))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# new merMod
Effect.merMod <- function(focal.predictors, mod, ..., KR=FALSE){
  if (KR && !requireNamespace("pbkrtest", quietly=TRUE)){
    KR <- FALSE
    warning("pbkrtest is not available, KR set to FALSE")}
  fam <- family(mod)
  args <- list(
    call = mod@call,
    coefficients = lme4::fixef(mod),
    family=fam,
    vcov = if (fam == "gaussian" && fam$link == "identity" && KR)
      as.matrix(pbkrtest::vcovAdj(mod)) else as.matrix(vcov(mod)))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# rlmer in robustlmm package
Effect.rlmerMod <- function(focal.predictors, mod, ...){
  args <- list(
    coefficients = lme4::fixef(mod),
    family=family(mod))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# clm in the ordinal package
Effect.clm <- function(focal.predictors, mod, ...){
  if (requireNamespace("MASS", quietly=TRUE)){
    polr <- MASS::polr}
  if(mod$link != "logit") stop("Effects only supports the logit link")
  if(mod$threshold != "flexible") stop("Effects only supports the flexible threshold")
  if(is.null(mod$Hessian)){
    message("\nRe-fitting to get Hessian\n")
    mod <- update(mod, Hess=TRUE)}
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  args <- list(
    type = "polr",
    coefficients = mod$beta,
    vcov = as.matrix(vcov(mod)[or, or]))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# clm2
Effect.clm2 <- function(focal.predictors, mod, ...){
  if (requireNamespace("MASS", quietly=TRUE)){
      polr <- MASS::polr}
  if(is.null(mod$Hessian)){
     message("\nRe-fitting to get Hessian\n")
     mod <- update(mod, Hess=TRUE)}
  if(mod$link != "logistic") stop("Effects only supports the logit link")
  if(mod$threshold != "flexible") stop("Effects only supports the flexible threshold")
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  args <- list(
    type = "polr",
    formula = mod$call$location,
    coefficients = mod$beta,
    vcov = as.matrix(vcov(mod)[or, or]))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

#clmm in ordinal package
Effect.clmm <- function(focal.predictors, mod, ...){
  if (requireNamespace("MASS", quietly=TRUE)){
    polr <- MASS::polr}
  if(is.null(mod$Hessian)){
    message("\nRe-fitting to get Hessian\n")
    mod <- update(mod, Hess=TRUE)}
  if(mod$link != "logit") stop("Only the logistic link is supported by Effects")
  if(mod$threshold != "flexible") stop("Only threshold='flexible supported by Effects")
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  skip <- length(unique(model.frame(mod)[,1])) - 1
  vcov <- matrix(NA, nrow=numBeta + skip, ncol=numBeta + skip)
  sel <- rownames(vcov(mod)) %in% names(mod$beta)
  vcov[1:numBeta, 1:numBeta] <- vcov(mod)[sel, sel]
  args <- list(
    type = "polr",
    formula = mod$formula,
    coefficients = mod$beta,
    vcov = as.matrix(vcov))
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
  args <- list(
    call = mod$call,
    formula = formula(mod),
    family=fam,
    coefficients = coef,
    vcov = vco)
  Effect.default(focal.predictors, mod, ..., sources=args)
}

