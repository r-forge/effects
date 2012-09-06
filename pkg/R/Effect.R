# Effect generic and methods
# John Fox and Sanford Weisberg
# last modified 2012-09-06 by J. Fox

Effect <- function(focal.predictors, mod, ...){
	UseMethod("Effect", mod)
}

Effect.lm <- function (focal.predictors, mod, xlevels = list(), default.levels = 10, given.values, 
		se = TRUE, confidence.level = 0.95, 
		transformation = list(link = family(mod)$linkfun, inverse = family(mod)$linkinv), 
		typical = mean, offset = mean, ...){
	if (missing(given.values)) 
		given.values <- NULL
	else if (!all(which <- names(given.values) %in% names(coef(mod)))) 
		stop("given.values (", names(given.values[!which]), ") not in the model")
	off <- if (is.numeric(offset) && length(offset) == 1) offset
			else if (is.function(offset)) {
				mod.off <- model.offset(model.frame(mod))
				if (is.null(mod.off)) 0 else offset(mod.off)
			}
			else stop("offset must be a function or a number")
	formula.rhs <- formula(mod)[[3]]
	model.components <- Analyze.model(focal.predictors, mod, xlevels, default.levels, formula.rhs)
	excluded.predictors <- model.components$excluded.predictors
	predict.data <- model.components$predict.data
	factor.levels <- model.components$factor.levels
	factor.cols <- model.components$factor.cols
	n.focal <- model.components$n.focal
	x <- model.components$x
	X.mod <- model.components$X.mod
	cnames <- model.components$cnames
	X <- model.components$X
	formula.rhs <- formula(mod)[c(1, 3)]
	Terms <- delete.response(terms(mod))
	mf <- model.frame(Terms, predict.data, xlev = factor.levels)
	mod.matrix <- model.matrix(formula.rhs, data = mf, contrasts.arg = mod$contrasts)
	wts <- mod$weights
	if (is.null(wts)) 
		wts <- rep(1, length(residuals(mod)))
	mod.matrix <- Fixup.model.matrix(mod, mod.matrix, model.matrix(mod), 
			X.mod, factor.cols, cnames, focal.predictors, excluded.predictors, typical, given.values)
	effect <- off + mod.matrix %*% mod$coefficients
	result <- list(term = paste(focal.predictors, collapse="*"), 
			formula = formula(mod), response = response.name(mod), 
			variables = x, fit = effect, x = predict.data[, 1:n.focal, drop=FALSE], model.matrix = mod.matrix, data = X, 
			discrepancy = 0, offset=off)
	if (se) {
		if (any(family(mod)$family == c("binomial", "poisson"))) {
			dispersion <- 1
			z <- qnorm(1 - (1 - confidence.level)/2)
		}
		else {
			dispersion <- sum(wts * mod$residuals^2)/mod$df.residual
			z <- qt(1 - (1 - confidence.level)/2, df = mod$df.residual)
		}
		V2 <- dispersion * summary.lm(mod)$cov
		V1 <- vcov(mod)
		V <- if (inherits(mod, "fakeglm")) 
					V1
				else V2
		vcov <- mod.matrix %*% V %*% t(mod.matrix)
		rownames(vcov) <- colnames(vcov) <- NULL
		var <- diag(vcov)
		result$vcov <- vcov
		result$se <- sqrt(var)
		result$lower <- effect - z * result$se
		result$upper <- effect + z * result$se
		result$confidence.level <- confidence.level
	}
	if (is.null(transformation$link) && is.null(transformation$inverse)) {
		transformation$link <- I
		transformation$inverse <- I
	}
	result$transformation <- transformation
	class(result) <- "eff"
	result
}

Effect.mer <- function(focal.predictors, mod, ...) {
    if (!require(lme4)) stop("the lme4 package is not installed")
    result <- Effect(focal.predictors, mer.to.glm(mod), ...)
    result$formula <- as.formula(formula(mod))
    result
}

Effect.lme <- function(focal.predictors, mod, ...) {
    if (!require(nlme)) stop("the nlme package is not installed")
    result <- Effect(focal.predictors, lme.to.glm(mod), ...)
    result$formula <- as.formula(formula(mod))
    result
}

Effect.gls <- function (focal.predictors, mod, xlevels = list(), default.levels = 10, given.values, 
                        se = TRUE, confidence.level = 0.95, 
                        transformation = NULL, 
                        typical = mean, ...){
    if (missing(given.values)) 
        given.values <- NULL
    else if (!all(which <- names(given.values) %in% names(coef(mod)))) 
        stop("given.values (", names(given.values[!which]), ") not in the model")
    formula.rhs <- formula(mod)[[3]]
    mod.lm <- lm(as.formula(mod$call$model), data=eval(mod$call$data))
    model.components <- Analyze.model(focal.predictors, mod.lm, xlevels, default.levels, formula.rhs)
    excluded.predictors <- model.components$excluded.predictors
    predict.data <- model.components$predict.data
    factor.levels <- model.components$factor.levels
    factor.cols <- model.components$factor.cols
    n.focal <- model.components$n.focal
    x <- model.components$x
    X.mod <- model.components$X.mod
    cnames <- model.components$cnames
    X <- model.components$X
    formula.rhs <- formula(mod)[c(1, 3)]
    nrow.X <- nrow(X)
    mf <- model.frame(formula.rhs, data=rbind(X[,names(predict.data),drop=FALSE], predict.data), 
                      xlev=factor.levels)
    mod.matrix.all <- model.matrix(formula.rhs, data=mf, contrasts.arg=mod$contrasts)
    mod.matrix <- mod.matrix.all[-(1:nrow.X),]
    mod.matrix <- Fixup.model.matrix(mod.lm, mod.matrix, model.matrix(mod.lm), 
                                     X.mod, factor.cols, cnames, focal.predictors, excluded.predictors, typical, given.values)
    fit.1 <- na.omit(predict(mod))
    mod.2 <- lm.fit(mod.matrix.all[1:nrow.X,], fit.1)
    class(mod.2) <- "lm"
    assign(".y", na.omit(model.response.gls(mod)), envir=.GlobalEnv)
    assign(".X", na.omit(mod.matrix.all[1:nrow.X,]), envir=.GlobalEnv)
    mod.3 <- update(mod, .y ~ .X - 1)
    remove(".X", ".y", envir=.GlobalEnv)
    discrepancy <- 100*mean(abs(fitted(mod.2)- fit.1)/(1e-10 + mean(abs(fit.1))))
    if (discrepancy > 1e-3) warning(paste("There is a discrepancy of", round(discrepancy, 3),
                                          "percent \n     in the 'safe' predictions used to generate effect", paste(focal.predictors, collapse="*")))
    effect <- mod.matrix %*% mod$coefficients
    result <- list(term = paste(focal.predictors, collapse="*"), 
                   formula = formula(mod), response = response.name(mod), 
                   variables = x, fit = effect, x = predict.data[, 1:n.focal, drop=FALSE], model.matrix = mod.matrix, data = X, 
                   discrepancy = discrepancy, offset=0)
    if (se){
        df.res <- mod$dims[["N"]] - mod$dims[["p"]]
        z <- qt(1 - (1 - confidence.level)/2, df=df.res)
        mod.2$terms <- terms(mod)
        V <- vcov(mod.3)
        vcov <- mod.matrix %*% V %*% t(mod.matrix)
        rownames(vcov) <- colnames(vcov) <- NULL
        var <- diag(vcov)
        result$vcov <- vcov
        result$se <- sqrt(var)        
        result$lower <- effect - z*result$se
        result$upper <- effect + z*result$se
        result$confidence.level <- confidence.level
    }
    if (is.null(transformation$link) && is.null(transformation$inverse)){
        transformation$link <- I
        transformation$inverse <- I
    }
    result$transformation <- transformation
    class(result) <- "eff"
    result
}
