################################################################
# An R script containing functions for model interpretation,
# assessment, and evaluation
################################################################

# TODO: make work for binary and gaussian responses
# Idea: Cross Validation takes in a gfo object
CrossValidate <- function(gfo, k.folds = 5, granularity = .005,
                          starter = list(), pre.test = FALSE) {

  # Currently set up for classification, binomial family and mids data only!
  #
  # Args:
  #   k.folds: Number of cross-validation folds
  #   granularity: For binary classification, the step size from 0 to 1 in
  #                computing evaluation metrics like precision
  #   starter: a list containing elements that allows restarting the
  #            cross validation if at least one fold was completed.
  #            Such a list is created at each iteration and placed in the
  #            global environment. It is called "starter.pack"
  #   pre.test: Whether to run the full model for every training split,
  #             to check for model stability before running a lengthy
  #             cross-validation run using backward elimination
  #
  # Returns:
  #  A gfo.cxv S3 object containing relevant outputs from the cross validation
  response <- all.vars(gfo$call.formula)[1]
  data <- gfo$data$data

  threshold.seq = seq(0, 1, granularity)
  if (length(starter) == 0) {
    split <- sample(1:k.folds, nrow(data), replace = TRUE)
    model.list <- list()
    min.redo.fold <- 1 # i.e., do them all
    precision.list <- list()
    accuracy.list <- list()
    recall.list <- list()
    complete.cases <- rep(0, length(threshold.seq))
  } else {
    model.list <- starter$model.list
    # re-label splits that are not in the current model list
    min.redo.fold <- min(k.folds - 1, length(model.list) + 1)
    split <- RelabelSplits(starter$split, c(min.redo.fold:k.folds))
    precision.list <- starter$precision.list
    accuracy.list <- starter$accuracy.list
    recall.list <- starter$recall.list
    complete.cases <- starter$complete.cases
  }
  if (pre.test) {
  # Optional pre-test that run the model on the full set of vars for each split
    options(width = 200)
    for (i in min.redo.fold:k.folds) {
      if (i == 0) rm(starter.pack, envir = .GlobalEnv)
      cat("##### Pre-test ", i, " of ", k.folds, " ######\n")
      temp <- FitModel(gfo$call.formula, data = gfo$data[split != i, ],
                       family = gfo$family)
      summary(temp)
    }
    rm(temp)
    check <- readline(prompt = "Check the runs for -999s. Continue? (y/n): ")
    if (check %in% c("n", "N")) stop("Unstable sample runs. Try Rerunning")
  }

  for (i in min.redo.fold:k.folds) {
    fold.str <- paste("\n########", "Fold", i, "of", k.folds, "#########\n")
    print(fold.str)
    train.df <- data[split != i, ]
    if (class(gfo$data) == "mids") {
      train.df <- mice(train.df, m = gfo$m, maxit = gfo$data$iteration,
                       predictorMatrix = gfo$data$predictorMatrix)
    }
    fold.model <- SequentiallyBuildModel(gfo$call.formula, train.df,
                                         gfo$var.select.cutoff, gfo$family,
                                         type = gfo$var.select.type)
    model.list[[i]] <- fold.model
    print(fold.str)
    test.df <- data[split == i, ]
    if (class(gfo$data) == "mids") {
      test.df <- mice(test.df, m = gfo$m, maxit = gfo$data$iteration,
                      predictorMatrix = gfo$data$predictorMatrix)
    }
      pred <- predict(fold.model, test.df)

    # TODO: THIS IS THE PART THATS BINARY-based
    fold.results.df <- data.frame(threshold = threshold.seq)
    fold.results.df$true.negatives <- NA
    fold.results.df$false.negatives <- NA
    fold.results.df$true.positives <- NA
    fold.results.df$false.positives <- NA

    for (row in 1:nrow(fold.results.df)) {
      flag <- as.numeric(pred >= fold.results.df[row, "threshold"])
      actual <- test.df[, response]
      fold.results.df[row, "true.negatives"] <- sum(flag == 0 & actual == 0)
      fold.results.df[row, "false.negatives"] <- sum(flag == 0 & actual == 1)
      fold.results.df[row, "true.positives"] <- sum(flag == 1 & actual == 1)
      fold.results.df[row, "false.positives"] <- sum(flag == 1 & actual == 0)
    }
    # Binary case accuracy, precision, recall
    fold.results.df <- within(fold.results.df, {
      accuracy <- (true.positives + true.negatives) /
                  (true.negatives + false.negatives + true.positives +
                   false.positives)
      precision <- true.positives / (true.positives + false.positives)
      recall <- true.positives / (true.positives + false.negatives)
    })
    accuracy.list[[i]] <- matrix(fold.results.df$accuracy, ncol = 1)
    recall.list[[i]] <- matrix(fold.results.df$recall, ncol = 1)
    precision.list[[i]] <- matrix(ifelse(is.nan(fold.results.df$precision), 0,
                                         fold.results.df$precision), ncol = 1)
    complete.cases <- complete.cases + !is.nan(fold.results.df$precision)
    starter.pack <<- list(split          = split,
                          model.list     = model.list,
                          accuracy.list  = accuracy.list,
                          recall.list    = recall.list,
                          precision.list = precision.list,
                          complete.cases = complete.cases)
  }
  metrics.df <- data.frame(threshold = threshold.seq)
  metrics.df$accuracy <- Reduce("+", accuracy.list) / k.folds
  metrics.df$precision <- Reduce("+", precision.list) / complete.cases
  metrics.df$recall <- Reduce("+", recall.list) / k.folds

  gfo.cxv <- list(gfo = gfo, metrics.df = metrics.df, cv.models = model.list,
                  cv.split = split, k.folds = k.folds)
  class(gfo.cxv) <- "gfo.cxv"
  return(gfo.cxv)
}

GetVariableImpacts <- function(obj, lower.quantile = .05,
                               upper.quantile = .95) {
  # Computes model extrapolations when fixed effects are modulated
  #
  # Args:
  #  obj: either a gfo.cxv (cross validation) or gfo object
  #  lower.quantile: the starting quantile of each continuous predictor
  #  upper.quantile: the ending quantile of each continuous predictor
  #
  # Returns: A data frame of the predicted outcomes under the modulations
  impact.df <- data.frame()
  if (class(obj) == "gfo.cxv") {
    gfo <- obj[["gfo"]]
    k.folds <- obj$k.folds
  } else if (class(obj) == "gfo") {
    gfo <- obj
    k.folds <- NA
    count <- NA
  } else {
    stop("Imput must be a gfo.cxv or gfo object")
  }

  if (class(gfo$data)[1] == "mids") {
    data <- gfo$data$data
  } else {
    data <- gfo$data
  }
  # TODO: Add functionality for ns() and cut() terms
  terms <- gfo$final.fixed.terms
  function.terms <- terms[grepl("\\(", terms)]
  for (term in terms) {
    if (class(obj) == "gfo.cxv") {
      count <- 0
      for (i in 1:k.folds) {
        final.cv.attribs <- attributes(terms(obj$cv.models[[i]]$formula))
        final.cv.terms <- final.cv.attribs$term.labels
        count <- count + as.numeric(term %in% final.cv.terms)
      }
    }
    term.sd <- gfo$sd.vars[term]
    if (term %in% function.terms) {
      term.obj <- data[, GetFnArg(term)]
    } else {
      term.obj <- data[, term]
    }
    cat("Getting impact for variable ", term, ", type: ", class(term.obj), "\n")
    if (is.numeric(term.obj) && !(term %in% function.terms)) {
      term.95 <- round(quantile(data[, term], upper.quantile, na.rm = TRUE), 3)
      term.05 <- round(quantile(data[, term], lower.quantile, na.rm = TRUE), 3)
      copy.df.95 <- data
      copy.df.05 <- data
      copy.df.95[, term] <- term.95
      copy.df.05[, term] <- term.05
      # Full model used to predict
      pred.high <- predict(gfo, copy.df.95)
      pred.low <- predict(gfo, copy.df.05)

      avg.diff <- round(mean(pred.high - pred.low, na.rm = TRUE), 3)
      avg.low <- round(mean(pred.low, na.rm = TRUE), 3)

    } else if (is.factor(term.obj)) {
      levels <- levels(term.obj)
      n.levels <- length(levels)
      cat(term, " is a factor with the following levels:\n")
      for (t in 1:n.levels) {
        cat("choice ", t, ": ", levels[t], "\n")
      }
      base.index <- readline(prompt = "Choose base level: ")
      cat("you chose ", base.index, "\n")
      base <- levels[as.numeric(base.index)]
      other.labels <- setdiff(levels, base)
      cat("base: ", base, "other.labels: ", other.labels, "\n")
      term.95 <- other.labels
      term.05 <- rep(base, n.levels - 1)
      pred.high <- c()
      pred.low <- c()
      avg.diff <- c()
      avg.low <- c()
      iter.count <- rep(count, n.levels - 1)
      for (k in 1:length(other.labels)) {
        print(term.95[k])
        print(term.05[k])
        copy.df.95 <- data
        copy.df.05 <- data
        copy.df.95[, term] <- factor(term.95[k], levels = levels)
        copy.df.05[, term] <- factor(term.05[k], levels = levels)
        pred.high <- c(pred.high, predict(gfo, copy.df.95))
        pred.low <- c(pred.low, predict(gfo, copy.df.05))
        avg.diff <- c(avg.diff, round(mean(pred.high - pred.low, na.rm = TRUE),
                                      3))
        avg.low <- c(avg.low, round(mean(pred.low, na.rm = TRUE), 3))
      }
    } else {
      print("Modulation of this term is currently not supported")
      term.05 <- NA
      term.95 <- NA
      pred.high <- NA
      pred.low <- NA
    }
    update <- data.frame(term = term, term.p.05 = as.character(term.05),
                         term.p.95 = as.character(term.95),
                         model.occurence.count = count,
                         model.occurance.pct = round(count / k.folds, 3),
                         avg.diff = avg.diff,
                         avg.low = avg.low)
    impact.df <- rbind(impact.df, update)
  }
  return(impact.df)
}

# idea: more manual, but you get what you want
GetVariableGroupImpact <- function(obj, low.val.list, high.val.list) {
  # Computes model extrapolations for a group of variables and prints output
  #
  # Args:
  #  obj: a gfo object
  #  low.val.list: a list with the names of the variables to be modulated and
  #    thier starting values
  #  high.val.list: a list with the names of the variables to be modulations
  #    and thier ending values
  #
  if (class(obj)[1] != "gfo") stop("Please supply a gfo object")
  if (class(obj$data)[1] == "data.frame") data <- obj$data
  if (class(obj$data)[1] == "mids") data <- complete(obj$data)

  if (!all(sort(names(low.val.list)) == sort(names(high.val.list)))) {
    stop("low.val.list and high value.list must have the same names")
  }
  var.vec <- names(low.val.list)

  copy.df.high <- data
  copy.df.low <- data
  for (term in var.vec) {
    copy.df.high[, term] <- high.val.list[[term]]
    copy.df.low[, term] <- low.val.list[[term]]
  }
  pred.high <- predict(obj, copy.df.high)
  pred.low <- predict(obj, copy.df.low)

  avg.diff <- round(mean(pred.high - pred.low, na.rm = TRUE), 3)
  avg.low <- round(mean(pred.low, na.rm = TRUE), 3)

  cat("-------------------------------------\n")
  cat(" Model Implied Group Impact Analysis \n")
  cat("    (Beware of extrapolation!)       \n")
  cat("-------------------------------------\n")

  cat('At the "low" values of\n')
  print(unlist(low.val.list))
  cat("the average predicted value is", avg.low, "\n\n")
  cat('The average delta from moving to the "high" values of\n')
  print(unlist(high.val.list))
  cat("is", avg.diff, "\n")

  return(list(pred.high=pred.high, pred.low=pred.low))
}
