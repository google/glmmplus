# Copyright 2013 Google Inc. All Rights Reserved.
# Author: Andy Chen <andych@google.com>
#         Ben Ogorek <baogorek@google.com>
# General prediction functions of lmer model.

predict.mer <- function(model, newdata = NULL, se.fit = FALSE,
                        n = 1000, verbose = FALSE, type = "response") {
  # Generate the variance of the predictors
  #
  # Args:
  #   model: a fitted lmer model object or a list containing the required
  #     lmer model elements.
  #   newdata: optionally, a data frame in which to look for variables with
  #     which to predict. If omitted, it returns the variance of the fitted
  #     lmer predictors.
  #   n: optionally, the number of iterations for simulating the predictors
  #     and get their variance.
  #   plot: optionally, a flag for generating the line plot for fits and
  #     variances.
  # Returns:
  #   The data frame with combined predicts and variance.
  # get fits
  if (is.null(newdata)) return(fitted(model))
  if (!is.data.frame(newdata)) stop("newdata should be a data.frame")

  thismodel <- ExtendMer(model, newdata)
  if (length(thismodel$ran.hat) == 0) ran.hat <- 0

  y.pred <- thismodel$invlinkfun(thismodel$fix.hat + thismodel$ran.hat)
  if (!se.fit) return(y.pred)

  y.sims <- SimulateMixedModel(model, newdata, n)

  obj <- list(fit   = y.pred, stdev = apply(y.sims, 1, sd),
              mean  = rowMeans(y.sims))
  class(obj) <- "predict.mer"
  return (obj)
}

plot.predict.mer <- function(y.pred, plot.smean = FALSE, plot.ylim = NULL,
                             xlab = NULL, ylab = NULL,
                             main = "lmer Predictions") {
  # Plot a line chart to show the predictions and their variances
  #
  # Args:
  #   y.pred: a data frame containing predicts and variances
  #   plot.smean: if TRUE, plot the mean of the simulation
  if (is.null(plot.ylim)) {
    y.range <- c(min(y.pred$fit - ifelse(is.na(y.pred$stdev), 0, y.pred$stdev)),
                 max(y.pred$fit + ifelse(is.na(y.pred$stdev), 0, y.pred$stdev)))
  } else {
    y.range <- plot.ylim
  }

  id <- 1:length(y.pred$fit)
  titlestr <- character(0)
  if (!is.null(xlab)) { titlestr <- paste(titlestr, ",xlab='", xlab, "'") }
  if (!is.null(ylab)) { titlestr <- paste(titlestr, ",ylab='", ylab, "'") }
  if (!is.null(main)) { titlestr <- paste(titlestr, ",main='", main, "'") }
  eval(parse(text = paste("plot(id, y.pred$fit, type='l', lty=1, col='blue',",
                          "ylim=y.range", titlestr, ")", sep="")))
  # plot(id, y.pred$fit, type='l', lty=1, col='blue', ylim=y.range,
  #      ylab='y.fit')
  if (plot.smean) {
    points(id, y.pred$mean, type='l', lty=2, col="darkgreen")
  }
  points(id, y.pred$fit - y.pred$stdev, type='l', lty=3, col="darkred")
  points(id, y.pred$fit + y.pred$stdev, type='l', lty=3, col="darkred")
  if (plot.smean) {
    legend('topright', c('y.fit', 'simu.mean', '+sigma', '-sigma'),
           col=c('blue', 'darkgreen', 'darkred', 'darkred'),
           lty = c(1, 2, 3, 3), pch = rep(21, 4))
  } else {
    legend('topright', c('y.fit', '+sigma', '-sigma'),
           col=c('blue', 'darkred', 'darkred'),
           lty = c(1, 3, 3), pch = rep(21, 3))
  }
}

SimulateMixedModel <- function(model, newdata, n = 1000) {
  # Generate the mean and sd of model predictions via simulation
  #
  # Args:
  #   model: a fitted lmer model object or a list containing the required
  #     lmer model elements.
  #   newdata: fresh data frame
  #   n: optionally, the number of iterations for simulating the predictors
  #     and get their variancel.
  #   plot: optionally, a flag for generating the line plot for fits and
  #     variances.
  # Returns:
  #   The data frame with combined predicts and variance.
  thismodel <- ExtendMer(model, newdata)
  lin.predictor <- thismodel$fix.hat + thismodel$ran.hat

  if (thismodel$family == "binomial") {
    y.sims <- simu.binom(n, 1000, plogis(lin.predictor)) / 1000
  } else if (thismodel$family == "poisson") {
    y.sims <- simu.pois(n, exp(lin.predictor))
  } else if (thismodel$family == "gaussian") {
    y.sims <- simu.norm(n, lin.predictor, thismodel$sigma.y)
  }
    y.pred <- data.frame(fit   = fitted(model),
                         stdev = apply(y.sims, 1, sd),
                         mean  = rowMeans(y.sims))
   return(y.pred)
}

ExtendMer <- function(model, newdata = NULL, verbose = FALSE) {
  # S3 (list) representation of an lmer object
  #
  # Args:
  #   model: a fitted lmer model object or a list containing the required
  #     lmer model elements.
  #   newdata: fresh data frame
  #   verbose: extra output
  # Returns:
  #   An S3 object of class "S3lmer" containing lmer object data
  thismodel <- list()
  thismodel$call <- model@call
  thismodel$terms <- terms(model)
  thismodel$ranef.list <- ranef(model)
  thismodel$fixef.list <- fixef(model)
  thismodel$fixef <- model@fixef
  thismodel$ranef <- model@ranef
  thismodel$varcorr <- VarCorr(model)
  class(thismodel) <- "S3lmer"

  if (length(grep("binomial", as.character(model@call))) == 1) {
    thismodel$family <- "binomial"
    thismodel$invlinkfun <- plogis
  } else if (length(grep("poisson", as.character(model@call))) == 1) {
    thismodel$family <- "poisson"
    thismodel$invlinkfun <- exp
  } else {
    thismodel$family <- "gaussian"
    thismodel$invlinkfun <- function(x) return(x)
  }

  if (is.null(newdata)) return(thismodel)

  modelZt <- get.modelZt(thismodel, newdata)  # ranef design matrix
  modelX <- get.modelX(thismodel, newdata)  # fixef design matrix
  thismodel$sigma.y <- get.lmer.sigma(model, 'data')

  thismodel$ran.ef <- numeric(dim(modelZt)[1])  # pre-allocate memory
  pos <- 0
  for (j in 1:length(ranef(model))) {  # dealing with all random effects
    sigma <- get.lmer.sigma(model, j)
    if (length(sigma) == 1) sigma <- sqrt(sigma)
    mu <- ranef(model)[[j]]
    iran.ef <- as.numeric(simu.norm(1, mu, sigma))
    thismodel$ran.ef[(pos + 1):(pos + length(iran.ef))] <- iran.ef
    pos <- pos + length(iran.ef)
  }

  thismodel$fix.hat <- modelX %*% model@fixef
  thismodel$ran.hat <- t(as.matrix(modelZt)) %*% thismodel$ranef

  return(thismodel)
}

get.modelZt <- function(model, newdata, verbose = FALSE) {
  # Get the model design matrix for the random efects
  #
  # Args:
  #   model: a fitted lmer model object or a list containing the required
  #     lmer model elements.
  #   newdata: optionally, a data frame in which to look for variables with
  #     which to predict.  If omitted, it returns the fitted lmer predictors.
  #   verbose: optionally, a flag for showing the dimensions of the model
  #     design matrix.
  # Returns:
  #   The design matrix for the random effects.
  thismodel <- model
  rn <- names(thismodel$ranef.list)

  # Set the starting length of each variable in the order of their positions
  # in the var list
  vlen <- numeric(length(rn))
  if (length(rn) > 1) {
    for (i in 2:length(rn)) {
      for (j in (1:(i - 1))) {
        d <- dim(thismodel$ranef.list[[j]])
        vlen[i] <- vlen[i - 1] + d[1] * d[2]
      }
    }
  }

  sparseMi <- numeric()
  sparseMj <- numeric()
  sparseX <- numeric()

  # For non-hierarchical random effects: assume all levels of the non-
  # hierarchical vars in the new data are in the training data; otherwise,
  # this section nees to be re-write for missing levels
  for (i in match(grep(':', rn, value=TRUE, invert=TRUE), rn)) {
    irn <- names(thismodel$ranef.list[[i]])
    for (j in 1:length(irn)) {  # for slopes
      isparseMi <- vlen[i] + (j - 1) * dim(thismodel$ranef.list[[i]])[1] +
        match(newdata[, rn[i]], rownames(thismodel$ranef.list[[i]]))
      # roll-up to the upper level for missing sub-level
      isparseMi <- isparseMi[!is.na(isparseMi)]
      if (length(isparseMi) > 0) {
        isparseMj <- 1:length(isparseMi)
      } else {
        isparseMj <- numeric()
      }

      # the grouping var rn[i] must be a factor
      if (irn[j] == '(Intercept)') {
        sparseX <- c(sparseX, rep(1, length(isparseMi)))
      } else if (is.factor(newdata[, irn[j]])) {
        sparseX <- c(sparseX, rep(1, length(isparseMi)))
      } else {
        sparseX <- c(sparseX, newdata[!is.na(isparseMi), irn[j]])
      }
      sparseMi <- c(sparseMi, isparseMi)
      sparseMj <- c(sparseMj, isparseMj)
    }
  }

  # for hierarchical random effects
  for (i in grep(':', rn)) {
    mixran.names <- gsub("[)]", "", gsub("[(]", "", split.str(rn[i], ":")))[[1]]
    nd <- newdata[, match(mixran.names, names(newdata))]

    cols <- dim(nd)[2]
    str <- character()
    sep = ":"
    for (j in 1:cols) { str <- paste(str, "nd[,", j, "],", sep="") }
    newrow <-
      eval(parse(text = paste("paste(", str, "sep='", sep, "')", sep="")))

    for (j in 1:length(irn)) {  # for slopes
      isparseMi <- vlen[i] + (j - 1) * dim(thismodel$ranef.list[[i]])[1] +
        match(newrow, rownames(thismodel$ranef.list[[i]]))
      if (sum(is.na(isparseMi)) > 0 && verbose) {
        print(paste("[", date(), "] warning: The levels of ", rn[i], ", ",
                    newrow[is.na(isparseMi)], ", are not in the model. ",
                    "Their effects were resplaced by their upper-level effects",
                    sep=""))
      }

      # roll-up to the upper level for missing sub-level
      isparseMi <- isparseMi[!is.na(isparseMi)]
      if (length(isparseMi) > 0) {
        isparseMj <- 1:length(isparseMi)
      } else {
        isparseMj <- numeric()
      }

      if (irn[j] == '(Intercept)') {  # the grouping var rn[i] must be a factor
        isparseX <- rep(1, length(isparseMi))
      } else if (is.factor(newdata[, irn[j]])) {
        isparseX <- rep(1, length(isparseMi))
      } else {
        isparseX <- newdata[which(!is.na(isparseMi)), irn[j]]
      }

      sparseX <- c(sparseX, isparseX)
      sparseMi <- c(sparseMi, isparseMi)
      sparseMj <- c(sparseMj, isparseMj)
    }
  }

  if (verbose) {
    print(paste("[", date(), "] messaging: ",
                "length(sparseMi)=", length(sparseMi), ", ",
                "length(sparseMj)=", length(sparseMj), ", ",
                "length(sparseX)=", length(sparseX), sep = ""))
  }
  ran.X <- Matrix::sparseMatrix(sparseMi, sparseMj, x = sparseX)

  return(ran.X)
}
get.modelX <- function(thismodel, newdata) {
  # Get the model design matrix for the fixed efects
  #
  # Args:
  #   thismodel: a fitted lmer model object or a list containing the required
  #     lmer model elements.
  #   newdata: optionally, a data frame in which to look for variables with
  #     which to predict.  If omitted, it returns the fitted lmer predictors.
  #
  # Returns:
  #   The design matrix for the fixed effects.
  fixed.effect.terms <- attr(thismodel$terms, "term.labels")
  fixed.form <- formula(paste("~", paste(fixed.effect.terms, collapse = "+")))
  X <- model.matrix(fixed.form, data = newdata)
  return(X)
}

simu.pois <- function(n, lambda) {
  # Simulate random numbers for possion distribution
  #
  # Args:
  #   n: number of observation groups.
  #   lambda: vector of (non-negative) means.
  # Returns:
  #   A matrix with the simulated results
  val <- matrix(0, nrow = length(lambda), ncol = n)
  for (i in 1:length(lambda)) {
    val[i, ] <- rpois(n, lambda[i])
  }
  return(val)
}

simu.binom <- function(n, size, prob) {
  # Simulate random numbers for binomial distribution
  #
  # Args:
  #   n: number of observation groups.
  #   size: number of trials.
  #   prob: probability vector of success on each trial.
  # Returns:
  #   A matrix with the simulated results

  val <- matrix(0, nrow = length(prob), ncol = n)
  for (i in 1:length(prob)) {
    val[i, ] <- rbinom(n, size, prob[i])
  }
  return(val)
}

simu.norm <- function(n, mu, sigma) {
  # Simulate random numbers for gaussian distribution
  #
  # Args:
  #   n: number of observation groups.
  #   mean: vector of means.
  #   sigma: vector of sigma. For multivariate vars, sigma is var-cor matrix.
  # Returns:
  #   A matrix with the simulated results

  if (length(sigma) > 1) {
    library(mvtnorm)
    if (nrow(mu) > 1 && ncol(mu) > 1) {
      # multiple row means using the same sigma
      val <- matrix(0, nrow = nrow(mu), ncol = n * ncol(sigma))
      for (icol in seq(1, n * ncol(sigma), by = ncol(sigma))) {
        ival <- matrix(0, nrow = nrow(mu), ncol = ncol(sigma))
        for (irow in 1:nrow(mu)) {
          ival[irow, ] <- rmvnorm(1, as.numeric(mu[irow, ]), sigma)
        }
        val[, icol:(icol + ncol(sigma) - 1)] <- ival
      }
    } else {
      val <- rmvnorm(n, as.numeric(mu), sigma)
    }
  } else { # for univariate, sigma is stdev
    val <- rnorm(n * nrow(mu), as.numeric(mu[, 1]), sigma)
    val <- matrix(val, nrow = nrow(mu), ncol = n)
  }

  colnames(val) <- rep(colnames(mu), n)
  return(val)
}

get.lmer.sigma <- function(model, var, stat = "sigma") {
  # Get the sigma matrix from the lmer model object.
  #
  # Args:
  #   model: lmer model object
  #   var: if it takes value "data", return the the scale parameter of from
  #        VarCorr; if takes value of an integer, retur ith component of the
  #        random effects from VarCorr function.
  #   stat: it takes value of "stdev" or "sigma", which returns the coresponding
  #         attributes of VarCorr function.
  # Returns:
  #   Sigma extracted from the lmer object.

  if (var == "data") {
    val <- attr(VarCorr(model), "sc")
    if (is.na(val)) val <- 1
    return(val)
  }

  if (is.numeric(var)) {
    id = var
  } else {
    varlist <- names(VarCorr(model))
    id <- which(varlist == var)
  }
  varcorr <- VarCorr(model)[[id]]

  if (stat == "stddev" || stat == "stdev") {
    val <- attr(varcorr, "stddev")
  } else if (substr(stat, 1, 4) == "corr") {
    val <- attr(varcorr, "correlation")
  } else if (stat == "sigma" || stat == "Sigma") {
    c <- attr(varcorr, "correlation")
    s <- attr(varcorr, "stddev")
    val <- c * s * rep(s, each = nrow(c))
  }

  val[is.na(val)] <- 0

  return(val)
}
