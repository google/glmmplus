# Copyright 2014 Google Inc. All rights reserved.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

coef2 <- function(object, ...) UseMethod('coef2', object)
coef2.glm <- function(object, ...) coef(object)
coef2.merMod <- function(object, ...) {
  # An S3 mer function resembling meant to resemble glm's coef
  #
  # Args:
  #  obj: a mer object
  #
  # Returns: A numeric vector with fixed effect coefficients
  val <- fixef(object)
  class(val) <- "coef.mer"
  return(val)
}

vcov2 <- function(object, ...) UseMethod('vcov2', object)
vcov2.glm <- function(object, ...) vcov(object)
vcov2.merMod <- function(object, ...) {
  # Overwriting mer variance function to match glm's format
  #
  # Args:
  #  obj: a mer object
  #
  #  Returns: The variance matrix as a native matrix object
#  rr <- as(lme4:::sigma(object)^2 *
#           chol2inv(object@RX, size = object@dims['p']), "dpoMatrix")
#  nms <- colnames(object@X)
#  dimnames(rr) <- list(nms, nms)
#  rr@factors$correlation <- as(rr, "corMatrix")
  rr <- vcov.merMod(object)
  var.mat <- try(as.matrix(rr))
  if (class(var.mat) == "matrix") {  # i.e., no errors
    rownames(var.mat) <- names(fixef(object))
    colnames(var.mat) <- names(fixef(object))
  }
  return(var.mat)
}

# Creating S3 generic function complete
complete <- function(df, ...) UseMethod('complete', df)
# compete is the identity function for data.frame inputs
complete.data.frame <- function(df, i = 1) df
# complete is the mice's complete function for mids objects
complete.mids <- function(df, i = 1) mice:::complete(df, i)

# The coef function for gfo objects is just a getter for the qbar element
coef.gfo <- function(obj) obj$qbar

scoef <- function(obj) {
  # Computes standardized coefficients of a gfo object
  # Args:
  #  obj: a gfo object
  #
  #  Returns: numerical vector of standardized coefficients
  response <- all.vars(obj$formula)[1]
  cont.preds <- intersect(names(obj$sd.vars), names(coef(obj)))
  denom <- 1
  if (obj$family()$family == "gaussian") {
    denom <- sd(complete(obj$data)[[response]])
  }
  result <- numeric(length(cont.preds))
  try(result <- coef(obj)[cont.preds] * obj$sd.vars[cont.preds] / denom)
  return(result)
}

vcov.gfo <- function(obj) {
  # Computes the Rubin's Rules-adjusted vaariance matrix of a gfo object
  # that provides row and column names
  # Args:
  #  obj: a gfo object
  #
  #  Returns: a dense variance matrix
  var.mat <- with(obj, ubar + (1 + 1 / m) * b)
  return(var.mat)
}

predict.gfo <- function(obj, newdata, type = "response") {
  # prediction function for a generalized fitted object
  #
  # Args:
  #  obj: a gfo object
  #  newdata: either a data.frame or a mids (multiply imputed data set) object
  #  type: TODO: this is not hooked up yet
  #
  #  Returns: numerical vector of fitted response values
  outer.pred.list <- list()
  test.m <- 1
  if (class(newdata) == "mids") test.m <- newdata$m
  for (i in 1:test.m) {
    inner.pred.list <- list()
    for (j in 1:obj$m) {
      model.j <- obj$fitted.list[[j]]
      inner.pred.list[[j]] <- as.matrix(predict(model.j, complete(newdata, i),
                                                type = "response"))
    }
    outer.pred.list[[i]] <- Reduce("+", inner.pred.list) / obj$m
  }
  pred <- Reduce("+", outer.pred.list) / test.m
  return(pred)
}

print.gfo <- function(obj) {
  # Printing function for a generalized fitted object (gfo)
  # Args:
  #  obj: a gfo object
  cat("Generalized Fitted Object 'gfo'\n\n")
  display.formula <- paste(deparse(obj$formula), collapse = " ")
  display.formula <- gsub("\\s+", " ", display.formula)
  cat("R Formula Syntax:", display.formula, "\n")
  if (class(obj$data) == "mids") cat("\nMultiple Imputation performed")
}

summary.gfo <- function(obj) {
  # Summary function for a generalized fitted object (gfo)
  # Args:
  #  obj: a gfo object
  print(obj)
  response <- all.vars(obj$formula)[1]
  summary.df <- data.frame("Intercept", NA,
                           round(coef(obj)["(Intercept)"], 3), NA)
  names(summary.df) <- c("Term", "Wald p-value", "Unstandardized coefs",
                         "Standarized coefs")
  if (identical(names(coef(obj)), "(Intercept)")) return("Intercept only model")
  for (term in names(obj$p.values)) {
    test.coefs <- obj$term.coef.map[[term]]
    p <- length(test.coefs)
    unstd.coefs <- round(coef(obj)[test.coefs], 3)
    std.coefs <- round(scoef(obj)[test.coefs], 3)
    pretty.p.val <- c(sprintf("%.4f", obj$p.values[term]), rep('pval above', p - 1))
    update <- data.frame(obj$term.coef.map[[term]], pretty.p.val, unstd.coefs, std.coefs)
    names(update) <- names(summary.df)
    summary.df <- rbind(summary.df, update)
  }
  if (obj$family()$family == "binomial") {
    summary.df$odds.ratio <- round(exp(summary.df[, "Unstandardized coefs"]), 3)
  }
  if (length(obj$random.terms) > 0) {
   cat("\n\n###  Random Effects  ###\n\n")
   var.comps <- list()
   for (i in 1:obj$m) {
     var.corr <- lme4:::formatVC(VarCorr(obj$fitted.list[[i]]))
     var.comps[[i]] <- as.numeric(var.corr[, "Variance"])
   }
   mean.vc <- Reduce("+", var.comps)
   var.corr[, 3] <- round(mean.vc, 2)
   var.corr[, 4] <- round(sqrt(mean.vc), 2)
   print(var.corr, quote = F, digits = 2) 
  }

  row.names(summary.df) <- NULL
  cat("\n### Type II analysis of fitted terms ###\n\n")
  print(summary.df)
  cat("\n ### Fit Statistics ###\n\n")
  # Idea: A batteries-included philosphy on output.
  if (obj$family()$family == "gaussian" && length(obj$random.terms) == 0) {
    cat("Well-defined R-squared is available\n")
    cat("R-squared:", round(obj$mcfaddens.r2, 4), "\n")
  } else if (obj$family()$family == "binomial") {
    cat("----Binomial response varible.----\n")
    cat("No well-defined R-squared is available.\n\n")
    cat("McFadden's generalized R-squared:", round(obj$mcfaddens.r2, 4), "\n\n")
    cat("Cox & Snell pseudo R-squared: ", round(obj$coxsnell.r2, 4), "\n")
  }
}

ComputeCoxSnellR2 <- function(fit, null.model, n) {
  # Computes the Cox & Snell R-squared
  #
  # Args:
  #  fit: a full model object, must support the LogLik function
  #  null.model: the reduced model object, must support the LogLike function
  #  n: the sample size
  #
  # Returns: the scalar Cox & Snell R-squared
  return(1 - (exp(logLik(null.model) - logLik(fit)))^(2 / n))
}

ComputeMcFaddensR2 <- function(fit, null.model) {
  # Computes the Cox & Snell R-squared
  #
  # Args:
  #  fit: a full model object, must support the LogLik function
  #  null.model: the reduced model object, must support the LogLike function
  #  n: the sample size
  #
  # Returns: the scalar Cox & Snell R-squared
  return(1 - logLik(fit) / logLik(null.model))
}

CollectModelQuantities <- function(fit, data, null.model) {
  # This is simply a device to avoid typing this multiple times
  obj <- list(fit.cov = vcov2(fit), fit.coef = coef2(fit),
              mcfaddens.r2 = ComputeCoxSnellR2(fit, null.model, nrow(data)),
              coxsnell.r2 = ComputeMcFaddensR2(fit, null.model),
              formula = formula, family = family)
  return(obj)
}

FitModel <- function(formula, data, family = gaussian) {
  # Fits model with, optionally, missing data and random effects
  #
  # Args:
  #  formula: a formula with optionally contains random effects
  #  data: a data.frame or mids object
  #  family: the family functions in R
  #
  # Returns: a gfo model object
  gfo <- BackwardEliminate(formula, data, 1, family, verbose = FALSE)
  gfo$var.select.type <- "none"
  gfo$var.select.cutoff <- NA
  return(gfo)
}

GetWaldPValue <- function(fitted, drop) {
  # Computes the Wald p-value based on the variance matrix of the estimates
  #
  # Args:
  #  fitted: a gfo object
  #  drop: a character vector of coefficient names that comprise the hypothesis
  #        of no effect (all zero coefficients)
  #
  # Return: scalar p-value
  if (norm(fitted$b) > 0) {  # uncertainty from missing data
    theta.bar <- fitted$qbar[drop]
    U.inv <- solve(fitted$ubar[drop, drop], tol = 1e-14)
    k <- length(theta.bar)
    r <- (1 + 1 / fitted$m) * tr(fitted$b[drop, drop] %*% U.inv) / k
    v <- k * (fitted$m - 1)
    w <- ifelse(v > 4, 4 + (v - 4) * (1 + (1 - 2 / v) / r)^2,
                0.5 * v * (1 + 1 / k) * (1 + 1 / r)^2)
    D <- (k * (1 + r))^(-1) * t(theta.bar) %*% U.inv %*% theta.bar
    p.value <- 1 - pf(D, k, w)
  } else {
    v <- vcov(fitted)[drop, drop]
    e <- matrix(coef(fitted)[drop], ncol = 1)
    df <- length(drop)
    d <- as.numeric(t(e) %*% solve(v) %*% e)
    p.value <- 1 - pchisq(d, df)
  }
  return(p.value)
}

GetEstimates <- function(data, ...) UseMethod("GetEstimates", data)

GetEstimates.data.frame <- function(data, formula, family, null.model,
                                    random.terms = character(0)) {
  # generalized fitted model (gfo) constructor for data.frame input data set
  #
  # Args:
  #  data: a data.frame
  #  formula: a formula that optionally contains random effects lme4 style
  #  family: one of the family functions (e.g., binomial), not in quotations
  #  null.model: the reduced model to be used in p-value computations
  #  random.terms: the vector of random terms from the model formula
  #
  #  Returns: a gfo S3 object
  if (length(random.terms) == 0) {
    model.fit <- glm(formula, family, data)
  } else {
    library(lme4)
    model.fit <- glmer(formula = formula, data = data, family = family)
  }
  quantities <- try(CollectModelQuantities(model.fit, data, null.model))
  analysis.list <- list()
  analysis.list[[1]] <- model.fit
  current.model <- list(ubar         = quantities[["fit.cov"]],
                        qbar         = quantities[["fit.coef"]],
                        b            = 0 * quantities[["fit.cov"]],
                        m            = 1,
                        mcfaddens.r2 = quantities[["mcfaddens.r2"]],
                        coxsnell.r2  = quantities[["coxsnell.r2"]],
                        formula      = formula, family = family,
                        fitted.list  = analysis.list,
                        sd.vars      = GetStandardDevs(data))
  class(current.model) <- "gfo"
  return(current.model)
}

GetEstimates.mids <- function(mids, formula, family, null.model,
                              random.terms = character(0)) {
  # generalized fitted model (gfo) constructor for mids (imputed data set)
  #
  # Args:
  #  data: a data.frame
  #  formula: a formula that optionally contains random effects lme4 style
  #  family: one of the family functions (e.g., binomial), not in quotations
  #  null.model: the reduced model to be used in p-value computations
  #  random.terms: the vector of random terms from the model formula
  #
  #  Returns: a gfo S3 object
  m <- mids$m
  library(multicore)
  if (length(random.terms) == 0) {
    analysis.list <- mclapply(c(1:m),
                              function(i) glm(formula, data = complete(mids, i),
                                              family = family),
                              mc.preschedule = TRUE, mc.cores = 5)
  } else {
    library(lme4)
    analysis.list <- mclapply(c(1:m),
                              function(i) glmer(formula,
                                               data = complete(mids, i),
                                               family = family),
                              mc.preschedule = FALSE, mc.cores = 5)
  }
  mat.list <- list()
  coef.list <- list()
  mcfaddens.list <- list()
  coxsnell.list <- list()
  l <- 0  # number of elements in the above list. Increments.
  for (i in 1:m) {
    current.analysis <- analysis.list[[i]]
    quantities <- try(CollectModelQuantities(current.analysis, complete(mids),
                                             null.model))
    if (class(quantities) == "try-error") {
      warning("Error. Model run dropped.")
      print(analysis.list[[i]])
    } else {
      l <- l + 1
      mat.list[[l]] <- quantities[["fit.cov"]]
      coef.list[[l]] <- quantities[["fit.coef"]]
      mcfaddens.list[[l]] <- quantities[["mcfaddens.r2"]]
      coxsnell.list[[l]] <- quantities[["coxsnell.r2"]]
    }
  }
  m <- length(mat.list)  # resetting m in case any analyses were dropped
  obj <- list(ubar = NA, qbar = NA, b = NA, m = m,
              mcfaddens.r2 = NA, coxsnell.r2 = NA,
              formula = formula, family = family, fitted.list = analysis.list,
              sd.vars = GetStandardDevs(mids$data))
  class(obj) = "gfo"
  if (m > 1) {
    obj[["qbar"]] <- Reduce("+", coef.list) / m
    obj[["ubar"]] <- Reduce("+", mat.list) / m
    obj[["mcfaddens.r2"]] <- mean(unlist(mcfaddens.list))
    obj[["coxsnell.r2"]] <- mean(unlist(coxsnell.list))
    b <- var(matrix(unlist(coef.list), nrow = m, byrow = TRUE))
    rownames(b) <- rownames(obj[["ubar"]])
    colnames(b) <- colnames(obj[["ubar"]])
    obj[["b"]] <- b
  }
  return(obj)
}
