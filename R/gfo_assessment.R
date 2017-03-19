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
############################################################

#' GetVariableImpacts
#'
#' Computes model extrapolations when fixed effects are modulated
#' 
#' @param obj a gfo object
#' @param lower.quantile the starting quantile of each continuous predictor
#' @param upper.quantile the ending quantile of each continuous predictor
#' 
#' @return A data frame of the predicted outcomes under the modulations
#'
#' @examples
#' my.mids <- ImputeData(testdata, m = 5, maxit = 5)
#' my.gfo <- BackwardEliminate(y ~ x + w + z, data = my.mids)
#' GetVariableImpacts(my.gfo)
#'
#' @export
GetVariableImpacts <- function(obj, lower.quantile = .05,
                               upper.quantile = .95) {
  impact.df <- data.frame()
  if (class(obj) == "gfo") {
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
  # TODO (baogorek): Add functionality for ns() and cut() terms
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

#' GetVariableGroupImpact: Model Implied Group Impact Analysis
#'
#' Shows how the predicted values of a model change as a set of variables
#' simultaneously change in value.
#' 
#' @param obj A gfo object created by FitModel, BackwardEliminate, or
#'                Forward Select
#' 
#' @param low.val.list A list mapping each variable name (from var.vec) to its
#'                    "low" value
#' 
#' @param high.val.list A list mapping each variable name (from var.vec) to its
#'                     "high" value
#' 
#' @examples 
#' my.gfo <- FitModel(y ~ x + w + z, data = testdata)
#' group.impacts <- GetVariableGroupImpact(my.gfo,
#'                                         low.val.list = list(x = -1, w = -1),
#'                                         high.val.list = list(x = 1, w = 1))
#'
#' @export
GetVariableGroupImpact <- function(obj, low.val.list, high.val.list) {
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
