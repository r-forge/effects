# Replication Code for Fox and Weisberg, "Visualizing Fit and Lack of Fit 
#   in Complex Regression Models with Predictor Effect Plots and Partial Residuals

## ------------------------------------------------------------------------
library("effects")

## ------------------------------------------------------------------------
m1 <- lm(log(infantMortality) ~ group * log(ppgdp), data = UN,
    subset = rownames(UN) != "Equatorial Guinea")
summary(m1)

## ----fig3-a,include=TRUE,fig.width=6.5,fig.height=5.5,fig.show='hide'----
plot(predictorEffect("group", m1,
        transformation = list(link = log, inverse = exp),
        xlevels = list(ppgdp = 10 ^ (2 : 5))),
    lines = list(multiline = TRUE),
    axes = list(
        x = list(rotate = 45),
        y = list(lab = "Infant Mortality", 
                 ticks = list(at = 2 ^ (1 : 8)))
        ),
    confint = list(style = "auto")
    )

## ----fig3-b,include=TRUE,fig.width=6.5,fig.height=5.5,fig.show='hide'----
plot(predictorEffect("ppgdp", m1,
        transformation = list(link = log, inverse = exp)),
    lines = list(multiline = TRUE),
    axes = list(
        x = list(rotate = 45),
        y = list(lab = "Infant Mortality", 
                 ticks = list(at = 2 ^ (1 : 8)))
        ),
    confint = list(style = "auto")
    )

## ----fig3a,include=TRUE,fig.width=8,fig.height=4,fig.show='hide'---------
m2 <- lm(infantMortality ~ group * ppgdp, data = UN)
plot(predictorEffects(m2, ~ ppgdp, partial.residuals = TRUE),
     axes = list(x = list(rotate = 25), y = list(lim = c(0, 150))),
     id = list(n = 1))

## ----fig3b,include=TRUE,fig.width=8,fig.height=4,fig.show='hide'---------
plot(predictorEffects(m1, ~ ppgdp, partial.residuals = TRUE),
    axes = list(x = list(rotate = 25)))

## ------------------------------------------------------------------------
summary(Cowles)

## ----include=FALSE-------------------------------------------------------
library("car") # for Anova(); we're suppressing messages that won't occur when the car package is updated to remove all data sets

## ----eval=FALSE----------------------------------------------------------
## library("car")

## ------------------------------------------------------------------------
mod.cowles.1 <- glm(volunteer ~ sex + neuroticism * extraversion,
    data = Cowles, family = binomial)
summary(mod.cowles.1)
Anova(mod.cowles.1)

## ----new-cowles-1,include=TRUE,fig.width=12,fig.height=4,fig.show='hide'----
plot(predictorEffects(mod.cowles.1,
        xlevels = list(extraversion = seq(0, 24, by = 6),
                     neuroticism = seq(0, 24, by = 6))),
    axes = list(y = list(type = "response")),
    lines = list(multiline = TRUE),
    rows = 1, cols = 3)

## ----new-cowles-2,include=TRUE,fig.width=12,fig.height=4,fig.show='hide'----
plot(predictorEffects(mod.cowles.1,
        ~ neuroticism, partial.residuals = TRUE),
     lattice = list(layout = c(4, 1)))

## ------------------------------------------------------------------------
library("splines")
mod.cowles.2 <- glm(volunteer ~
        sex + ns(neuroticism, 5) * ns(extraversion, 5),
    data = Cowles, family = binomial)
anova(mod.cowles.1, mod.cowles.2, test = "Chisq")
cbind(AIC(mod.cowles.1, mod.cowles.2),
      BIC(mod.cowles.1, mod.cowles.2))

## ------------------------------------------------------------------------
summary(Prestige)

## ------------------------------------------------------------------------
Prestige$type <- factor(Prestige$type,
    levels = c("bc", "wc", "prof"))
mod.prestige.1 <- lm(prestige ~ income + education + type,
    data = Prestige)
summary(mod.prestige.1)
Anova(mod.prestige.1)

## ----fig-prestige-1,include=TRUE,fig.width=5,fig.height=5,fig.show='hide'----
plot(predictorEffects(mod.prestige.1, ~ income,
    partial.residuals = TRUE))

## ----fig-prestige-2,include=TRUE,fig.width=8,fig.height=4,fig.show='hide'----
plot(Effect(c("income", "type"), mod.prestige.1,
        partial.residuals = TRUE),
    partial.residuals = list(span = 0.9),
    axes = list(x = list(rotate = 25)),
    lattice = list(layout = c(3, 1)))

## ------------------------------------------------------------------------
mod.prestige.2 <- lm(prestige ~ type * income + education,
    data = Prestige)
anova(mod.prestige.1, mod.prestige.2)

## ----fig-prestige-3,include=TRUE,fig.width=8,fig.height=4,fig.show='hide'----
plot(Effect(c("income", "type"), mod.prestige.2,
        partial.residuals = TRUE),
    partial.residuals = list(span = 0.9),
    axes = list(x = list(rotate = 25)),
    lattice = list(layout = c(3, 1)))

## ----include=FALSE-------------------------------------------------------
mvrunif <- function(n, R, min = 0, max = 1){
# method (but not code) from E. Schumann,
# "Generating Correlated Uniform Variates"
# URL:
# <http://comisef.wikidot.com/tutorial:correlateduniformvariates>
# downloaded 2015-05-21
if (!is.matrix(R) || nrow(R) != ncol(R) ||
max(abs(R - t(R))) > sqrt(.Machine$double.eps))
stop("R must be a square symmetric matrix")
if (any(eigen(R, only.values = TRUE)$values <= 0))
stop("R must be positive-definite")
if (any(abs(R) - 1 > sqrt(.Machine$double.eps)))
stop("R must be a correlation matrix")
m <- nrow(R)
R <- 2 * sin(pi * R / 6)
X <- matrix(rnorm(n * m), n, m)
X <- X %*% chol(R)
X <- pnorm(X)
min + X * (max - min)
}

gendata <- function(n = 5000, R, min = -2, max = 2, s = 1.5,
model = expression(x1 + x2 + x3)){
data <- mvrunif(n = n, min = min, max = max, R = R)
colnames(data) <- c("x1", "x2", "x3")
data <- as.data.frame(data)
data$error <- s * rnorm(n)
data$y <- with(data, eval(model) + error)
data
}

R <- function(offdiag = 0, m = 3){
R <- diag(1, m)
R[lower.tri(R)] <- R[upper.tri(R)] <- offdiag
R
}

## ----include=FALSE-------------------------------------------------------
set.seed(682626)
Data.4 <- gendata(R = R(0.5), model = expression(x1 ^ 2 + x2 * x3))
mod.4 <- lm(y ~ x1 + x2 + x3, data = Data.4)

## ----fig-contrived-4a,include=FALSE,fig.width=12,fig.height=4,fig.show='hide'----
plot(predictorEffects(mod.4, partial.residuals=TRUE),
     partial.residual = list(pch = ".", col = "#FF00FF80"),
     axes = list(x = list(rotate = 45)),
     rows = 1, cols = 3)

## ----fig-contrived-4b,include=FALSE,fig.width=12,fig.height=4,fig.show='hide'----
plot(Effect(c("x2", "x3"), mod.4, partial.residuals = TRUE),
     partial.residual = list(pch = ".", col = "#FF00FF80"),
      axes = list(x = list(rotate = 45)),
     lattice = list(layout = c(4, 1)))

## ----fig-contrived-4c,include=FALSE,fig.width=12,fig.height=4,fig.show='hide'----
plot(Effect(c("x1", "x2"), mod.4, partial.residuals = TRUE),
    partial.residual = list(pch = ".", col = "#FF00FF80"),
    axes = list(x = list(rotate = 45)),
    lattice = list(layout = c(4, 1)))

## ----fig-contrived-5a,include=FALSE,fig.width=5,fig.height=4,fig.show='hide'----
mod.5 <- lm(y ~ poly(x1, 2) + x2 * x3, data = Data.4)
plot(Effect("x1", mod.5, partial.residuals = TRUE),
     partial.residual = list(pch = ".", col = "#FF00FF80", span = 0.2))

## ----fig-contrived-5b,include=FALSE,fig.width=12,fig.height=4,fig.show='hide'----
plot(Effect(c("x2", "x3"), mod.5, partial.residuals = TRUE),
     partial.residual = list(pch = ".", col = "#FF00FF80"),
     axes = list(x = list(rotate = 45)),
     lattice = list(layout = c(4, 1)), span = 0.5)

## ----fig-contrived-5c,include=FALSE,fig.width=12,fig.height=4,fig.show='hide'----
plot(Effect(c("x1", "x2"), mod.5, partial.residuals = TRUE),
    partial.residual = list(pch = ".", col = "#FF00FF80", span = 0.35),
    axes = list(x = list(rotate = 45)),
    lattice = list(layout = c(4, 1)))

