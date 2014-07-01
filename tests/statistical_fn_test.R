#library(mice)
#library(lme4)
#library(RUnit)
#library(splines)
#source("../R/analysis_functions.R");
#source("../R/helper_functions.R");
#source("../R/sequential.R")
#source("../R/lmer_prediction.R")
#source("../R/gfo_assessment.R")
#source("../R/mice_extentions.R")

set.seed(1234)
N <- 500
complete.df <- data.frame(x = rnorm(N), z = rnorm(N), w = rnorm(N),
                          factor.1 = factor(sample(30, N, replace = TRUE)),
                          factor.2 = factor(sample(5, N, replace = TRUE)))
random.effects.1 <- rnorm(30, mean = 0, sd = 2)
random.effects.2 <- rnorm(5, mean = 0, sd = 1.5)
# factor.2 used as both a random effect and a fixed categorical factor
complete.df$y <- with(complete.df, random.effects.1[factor.1] +
                 random.effects.2[factor.2] + 0.9 * sin(x)
     + 1.0 * w + 2 * rnorm(N))

complete.df$y.binary <- as.numeric(complete.df$y > 0)

spline.fn <- MakeSplineStr("x", complete.df, interior.knots = c(-.5, .5))

missing.df <- complete.df
missing.df[sample(1:N, round(N / 4)), "x"] <- NA
missing.df[sample(1:N, round(N / 4)), "z"] <- NA

test.mids <- ImputeData(missing.df, m = 3, maxit = 3,
                        droplist = c("y.binary", "factor.1"))
rm(complete.df)
# Creating package data set
# testdata <- missing.df
# save(list = "testdata", file = "../data/testdata.rda", compress = TRUE)

TestWaldPValues <- function() {
  # Complete case
  test.glm <- glm(y ~ x, data = missing.df)
  glm.p.value <- summary(test.glm)$coefficients["x", "Pr(>|t|)"]

  # Manually create a gfo object below
  test.gfo <- list(qbar = coef(test.glm), ubar = vcov(test.glm),
                   b = 0 * vcov(test.glm), m = 1)
  class(test.gfo) <- "gfo"

  # The wald p-value shouldn't be too different glm's default
  glmmplus.pvalue <- glmmplus:::GetWaldPValue(test.gfo, "x")
  checkTrue(abs(glm.p.value - glmmplus.pvalue) < .01)
  null.model <- glm(y ~ 1, data = missing.df)

  # mira: "multiply imputed repeated analysis"
  test.mira.full <- with(test.mids, glm(y ~ w + z + x))
  test.mira.red <- with(test.mids, glm(y ~ w + z))
  mice.p.value <- pool.compare(test.mira.full, test.mira.red)$pvalue

  # Our alternative gfo object
  # Note that interally to glmmplus, GetEstimates() is generic
  glmmplus.pooled <- glmmplus:::GetEstimates.mids(test.mids, y ~ w + z + x,
                                                  gaussian, null.model)
  test.coef.names <- glmmplus:::GetCoefNames("x", complete(test.mids))

  # How did we do matching mice on a nested hypothesis p-value?
  glmmplus.p.value <- glmmplus:::GetWaldPValue(glmmplus.pooled, test.coef.names)
  checkEqualsNumeric(log(mice.p.value), log(glmmplus.p.value))

  # How did we do versus mice's pooled analysis?
  test.mipo <- pool(test.mira.full) # mice's pooled analysis
  checkEqualsNumeric(glmmplus.pooled$ubar, test.mipo$ubar)
  checkEqualsNumeric(glmmplus.pooled$qbar, test.mipo$qbar)
  checkEqualsNumeric(glmmplus.pooled$b, test.mipo$b)
  rm(glmmplus.pooled)
  rm(test.mipo)
  rm(test.glm)
  rm(test.gfo)
  rm(test.mira.full)
  rm(test.mira.red)
}

TestSequentialEdgeCases <- function() {
  # Edge Case 1 - no terms remain
  glm.null <- BackwardEliminate(y ~ z, missing.df)
  checkEquals(length(glm.null$qbar[-1]), 0)

  glm.null <- ForwardSelect(y ~ z, missing.df)
  checkEquals(length(glm.null$qbar[-1]), 0)

  glm.null <- BackwardEliminate(y.binary ~ z, missing.df, family = binomial)
  checkEquals(length(glm.null$qbar[-1]), 0)
  rm(glm.null)
  # Edge Case 2 - all terms remain
  glm.complete <- BackwardEliminate(y ~ x + w, missing.df)
  checkEquals(sort(names(glm.complete$qbar[-1])), c("w", "x"))

  glm.complete <- ForwardSelect(y ~ x + w, missing.df)
  checkEquals(sort(names(glm.complete$qbar[-1])), c("w", "x"))
  rm(glm.complete)
  lmer.complete <- BackwardEliminate(y ~ (1 | factor.1) + x + w, missing.df)
  checkEquals(sort(names(lmer.complete$qbar[-1])), c("w", "x"))
  rm(lmer.complete)
}

TestImputationGLMM <- function() {
  # glm all linear terms, multiply imputed data set
  glm.impute <- BackwardEliminate(y ~ x + w + z, test.mids)
  checkEquals(sort(names(glm.impute$qbar[-1])), c("w", "x"))

  glm.impute <- BackwardEliminate(y.binary ~ x + w + z, test.mids,
                                  family = binomial)
  checkEquals(sort(names(glm.impute$qbar[-1])), c("w", "x"))
  rm(glm.impute)

  # lmer, linear terms, 2 random effects, complete data
  lmer.impute <- BackwardEliminate(y ~ (1 | factor.1) + (1 | factor.2) +
                                   x + w + z, test.mids)
  # Making sure the formula ended up correct
  checkEquals(lmer.impute$formula,
              formula("y ~ (1 | factor.1) + (1 | factor.2) + x + w"))
  checkEquals(sort(names(lmer.impute$qbar[-1])), c("w", "x"))
  rm(lmer.impute)
}

TestGLMMsWithBasisFunctions <- function() {
  # Instead of variable names, looking for first letter of ns() and cut()
  form <- formula(paste("y.binary ~ (1 | factor.1) +", spline.fn, "+ w + z"))
  lmer.spline <- ForwardSelect(form, test.mids, family = binomial)
  checkTrue(all(substr(names(lmer.spline$qbar[-1]), 1, 1) %in% c("n", "w")))
  rm(lmer.spline)

  lmer.cut <- BackwardEliminate(y.binary ~ (1 | factor.1) +
                                cut(x, c(-10, -1, 0, 1, 10)) + w + z,
                                test.mids, family = binomial)
  checkTrue(all(substr(names(lmer.cut$qbar[-1]), 1, 1) %in% c("c", "w")))
  rm(lmer.cut)
}