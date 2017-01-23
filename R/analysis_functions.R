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
  val <- fixef(object)
  class(val) <- "coef.mer"
  return(val)
}
coef2.lme <- function(object, ...) {
  return(coef2.merMod(object))
}

vcov2 <- function(object, ...) UseMethod('vcov2', object)
vcov2.glm <- function(object, ...) vcov(object)
vcov2.lme <- function(object, ...) vcov(object)
vcov2.merMod <- function(object, ...) {
  # Overwriting mer variance function to match glm's format
  #
  # Args:
  #  obj: a mer object
  #
  #  Returns: The variance matrix as a native matrix object
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
  if (class(obj$data) == "mids") {
    cat("\nNote: Multiple Imputation of missing values performed\n")
    cat("P-values adjusted according to Rubin's Rules\n\n")
  }
}

plotFSR <- function(gfo) {
  # Plot how the False Selection Rate changes during the selection techniqu
  # TODO (baogorek): need a stop() for when variable selection was not used
  term <- c()
  formula <- c()
  fsr <- c()

  for (i in 1:length(gfo$history)) {
    term <- c(term, gfo$history[[i]]$term)
    fsr <- c(fsr, gfo$history[[i]]$fsr)
  }
  par(las = 2) # make label text perpendicular to axis
  par(mar = c(5, 11, 4, 2)) #  margin: bottom, left, top, right 
  barplot(fsr, names.arg = substr(term,1, 27), horiz = TRUE,
          main = "Fast False Selection Rate Estimation",
          ylab = "", xlab = "Estimated FSR", xlim = c(0, .5))
  # TODO (baogorek): the 12 in arrows and 6 in text need to be dynamic
  arrows(.4, 1, .4, 12)
  text(.4, 6, "iteration", pos = 4)
  if (gfo$var.select.type == "forward") {
    message <- "term added from forward selection"
  } else if (gfo$var.select.type == "backward") {
    message <- "term dropped from backward elimination"
  }
  par(las = 3)
  # TODO (baogorek): make padj dynamic (if necessary)
  mtext(message, 2, col = "black", padj = -13.5, cex = 1.2, font = 3)
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
    update <- data.frame(obj$term.coef.map[[term]], pretty.p.val, unstd.coefs,
                         std.coefs)
    names(update) <- names(summary.df)
    summary.df <- rbind(summary.df, update)
  }
  if (obj$family()$family == "binomial") {
    summary.df$odds.ratio <- round(exp(summary.df[, "Unstandardized coefs"]), 3)
  }
  if (length(obj$random.terms) > 0) {
    var.comps <- list()
    ar.vec <- c()
    nlme.flag <- FALSE
    for (i in 1:obj$m) {
      if (class(obj$fitted[[i]]) == "lme") {
        nlme.flag <- TRUE
        library(nlme)
        var.comps[[i]] <- c(getVarCov(obj$fitted.list[[i]]),
                            obj$fitted.list[[i]]$sigma^2)
        # TODO un-hard code subject
        var.corr <- data.frame(Groups = c("Subject",  "Residual"),
                               Name = c("(Intercept)", ""),
                               Variance = numeric(2))
        model.struct <- summary(obj$fitted[[i]])$modelStruct
        cor.struct <- summary(model.struct)$corStruct
        # Changing class to corAR1
        #class(cor.struct) <- attr(cor.struct, "oClass")
        #print(as.numeric(coef(cor.struct, unconstrained = FALSE)))
        ar.parm <- nlme:::coef.corAR1(cor.struct, unconstrained = FALSE)
        ar.vec <- c(ar.vec, ar.parm)
      } else {
        # var.corr contains names that do not change with i
        var.corr <- lme4:::formatVC(VarCorr(obj$fitted.list[[i]]))
        var.comps[[i]] <- (as.numeric(var.corr[, "Std.Dev."]))^2
      }
    }
    mean.vc <- Reduce("+", var.comps) / length(var.comps)
    var.corr[, 3] <- round(sqrt(mean.vc), 2)
    cat("\n\n###  Random Effects  ###\n\n")
    print(var.corr, quote = F, digits = 2) 
    if (nlme.flag) {
      cat("\n AR1 parameter estimate from nlme:\n")
      cat(round(mean(ar.vec), 2), "\n")
    } 
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
  } else {
    cat("No fit statistics have been selected for this model\n")
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

FitModel <- function(formula, data, family = gaussian, ts.model = NULL) {
  # Fits model with, optionally, missing data and random effects
  #
  # Args:
  #  formula: a formula with optionally contains random effects
  #  data: a data.frame or mids object
  #  family: the family functions in R
  #
  # Returns: a gfo model object
  gfo <- BackwardEliminate(formula, data, 1, family, ts.model, verbose = FALSE)
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
                                    random.terms = character(0),
                                    ts.model = NULL) {
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
    library(splines)
    model.fit <- glm(formula = formula, family = family(), data = data)
    if (!is.null(ts.model)) {
      stop("Implement subject-specific random effect to use ts.model")
    }
  } else {
    library(lme4)
    library(splines)
    if (family()$family == "gaussian") {
      if (!is.null(ts.model)) {
        if (length(random.terms) > 1) {
          stop("Only one random term is allowed for use with ts.model")
        }
        library(nlme)
        # TODO: fix D.R.Y. violation (sequential.R)
        all.terms <- attributes(terms(formula))$term.labels
        fixed.terms <- all.terms[!grepl("\\|", all.terms)]
        response <- all.vars(formula)[1]
        fixed.formula <- as.formula(paste0(response, "~",
                                paste(fixed.terms, collapse = "+")))
        random.formula <- as.formula(paste0("~", random.terms[1]))
        if (ts.model == "ar1") {
          model.fit <- lme(fixed.formula, data = data,
                           random = random.formula,
                           cor = corAR1(0.5, form = random.formula))
        } else {
          stop("That ts.model option has not yet been implemented")
        }
      } else {
        model.fit <- lmer(formula = formula, data = data)
      }
    } else {
      if (!is.null(ts.model)) {
        stop("Sorry, within-subject time series only available for gaussian")
      }
      model.fit <- glmer(formula = formula, data = data, family = family)
    }
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
                              random.terms = character(0), ts.model = NULL) {
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
  library(parallel)
  if (length(random.terms) == 0) {
    library(splines)
    analysis.list <- lapply(c(1:m),
                              function(i){
                                library(splines)
                                glm(formula = formula, data = complete(mids, i),
                                    family = family)})
                     
  } else {
    if (family()$family == "gaussian") {
      # Handle Gaussian ts models
      if (!is.null(ts.model)) {
        if (length(random.terms) > 1) {
          stop("Only one random term is allowed for use with ts.model")
        }
        all.terms <- attributes(terms(formula))$term.labels
        fixed.terms <- all.terms[!grepl("\\|", all.terms)]
        response <- all.vars(formula)[1]
        fixed.formula <- as.formula(paste0(response, "~",
                                paste(fixed.terms, collapse = "+")))
        random.formula <- as.formula(paste0("~", random.terms[1])) 
        if (ts.model == "ar1") {
          library(nlme)
          library(splines)
          analysis.list <- lapply(c(1:m), function(i) {
                           library(nlme)
                           library(splines)
                           return(lme(fixed.formula, data = complete(mids, i),
                                      random = random.formula,
                                      cor = corAR1(0.5, form = random.formula)))
                                  }) 
        } else {
          stop("That ts.model option has not yet been implemented")
        }
      } else {
      # Handle Gaussian non-ts models 
        library(lme4)
        library(splines)
        analysis.list <- lapply(c(1:m),
                               function(i) {
                                 library(lme4)
                                 library(splines)
                                 return(lmer(formula = formula,
                                             data = complete(mids, i)))})
      }
    } else {
       # Handle non-Gaussian random effects models 
        if (!is.null(ts.model)) {
               stop("Sorry, within-subject time series only available for gaussian")
        }
        library(lme4)
        library(splines)
        analysis.list <- lapply(c(1:m),
                               function(i){
                                 library(lme4)
                                 library(splines)
                                 return(glmer(formula = formula,
                                        data = complete(mids, i),
                                        family = family))
                                })
     }
  }
  mat.list <- list()
  coef.list <- list()
  mcfaddens.list <- list()
  coxsnell.list <- list()
  l <- 0  # number of elements in the above list. Increments.
  for (i in 1:m) {
    current.analysis <- analysis.list[[i]]
    quantities <- try(CollectModelQuantities(current.analysis, complete(mids, i),
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
