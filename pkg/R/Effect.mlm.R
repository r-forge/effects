#' Calculate Effects for term(s) in a Multivariate Linear Model 

#' This function calculates effects for some or all of the response
#' variables for a term in a multivariate linear model.
#' It uses \code{\link[stats]{update}} to evaluate the linear model
#' for each response variable.

#' @param focal.predictors the name(s) of predictor variables in the \code{mod} for which
#'        to calculate effects
#' @param mod a \code{mlm} object
#' @param response name(s) or indices of one or more response variable(s).  The default is to use
#'        all responses
#' @return For one response, an \code{eff} object, otherwise a class \code{efflist} object, 
#'          containing one \code{eff} object for each \code{response}

#' @note This may require some minor adjustments to print and summary methods related to
#'       labeling of output (response and term)


Effect.mlm <- function(focal.predictors, mod, response, ...) {
	if (missing(response)) {
		mod.frame <- model.frame(mod)
    response <- colnames(model.response(mod.frame))
	}
	else if (is.numeric(response)) {
		mod.frame <- model.frame(mod)
    response.names <- colnames(model.response(mod.frame))
    response <- response.names[response]
	}
	
	if (length(response)==1) {
			mod.1 <- update(mod, as.formula(paste(response, " ~ .")))
			result <- Effect(focal.predictors, mod.1,  ...)
	}
	else {
		result <- as.list(NULL)
		for (resp in response) {
			mod.1 <- update(mod, as.formula(paste(resp, " ~ .")))
			lab <- resp
			result[[lab]] <- Effect(focal.predictors, mod.1,  ...)
		}
		class(result) <- "efflist"
	}
	result
}

