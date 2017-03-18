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

GetFastFSR <- function(n.total.vars, n.model.vars, alpha, verbose = FALSE) {
  # Compute Fast FSR estimate from
  # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2748177/
  #
  # Args:
  #  n.total.vars: Number of variables eligible for model
  #  n.model.vars: Number of variables kept in model
  #  alpha: the traditional alpha-to-enter or alpha-to-stay
  #
  # Returns: The FSR estimate of Boos, Stefanski, and Wu (2009)
  if (verbose) {
  cat ("Echoing Fast FSR paramters: \n n.total.vars: ", n.total.vars,
       "n.model.vars: ", n.model.vars, ", alpha: ", alpha, "\n")
  }
  N.hat <- n.total.vars - n.model.vars  # Estimated number of unimportant vars
  theta.hat <- alpha  # estimated rate of unimportant variables entering model
  return(N.hat * theta.hat / (1 + n.model.vars))
}

#' Backward Elimination for Generalized (Mixed) Linear Models with Missing Data
#'
#' This is a backwards elimination procedure for generic procedures that
#' ouput p-values
#' 
#' 
#' @param formula A formula which may contain random effects according to the
#'                lme4 package's specification.
#' @param data    Either a mids object from the mice package, or a data frame.
#' @param cutoff The alpha level which determines the stopping rule. Once all
#'               remaining model terms fall below this value,
#'               the procedure terminates.
#' @param family Any family accepted by glm or lmer. Do not use quotation marks.
#' 
#' @return A gfo object
#' 
#' @examples
#' # A sample data set with testdata values
#' head(testdata)
#' 
#' # creating a Muliply Imputed Data Set (mids) object
#' my.mids <- ImputeData(testdata, m = 5, maxit = 5)
#' 
#' # a single imputation
#' complete1 <- complete(my.mids)
#' 
#' # Backwards elimination for fixed effect models
#' BackwardEliminate(y ~ x + w + z, data = complete1)
#' BackwardEliminate(y ~ x + w + z, data = my.mids)
#' 
#' # Backwards elimination for mixed (fixed and random) models
#' BackwardEliminate(y ~ (1 | factor.1) + x + w + z, data = complete1)
#' BackwardEliminate(y ~ (1 | factor.1) + x + w + z, data = my.mids)
#'
#' @references
#' Douglas Bates and Martin Maechler (2010). lme4: Linear mixed-effects models
#' using S4 classes. R package version 0.999375-37. http://CRAN.R-project.org/package=lme4
#' 
#' Stef van Buuren, Karin Groothuis-Oudshoorn (2011).
#' mice: Multivariate Imputation by Chained Equations in R.
#'  Journal of Statistical Software, 45(3), 1-67. URL http://www.jstatsoft.org/v45/i03/.
#' 
#' @export
BackwardEliminate <- function(formula, data, cutoff = .05,
                              family = gaussian, ts.model = NULL,
                              verbose = TRUE) {
  final.model <- SequentiallyBuildModel(formula, data, cutoff, family,
                                        ts.model, type = "backward", verbose)
  return(final.model)
}

#' Forward Select
#'
#' Selection for Generalized (Mixed) Linear Models with Missing Data
#' 
#' @param formula A formula which may contain random effects according to the
#'                lme4 package's specification.
#' @param data Either a mids object from the mice package, or a data frame.
#' @param cutoff The alpha level which determines the stopping rule. Once all
#'               remaining model terms fall below this value,
#'               the procedure terminates.
#' @param family Any family accepted by glm or lmer. Do not use quotation
#'                 marks.
#'
#' @examples
#' # A sample data set with missing values
#' head(testdata)
#' 
#' # creating a Muliply Imputed Data Set (mids) object
#' my.mids <- ImputeData(testdata, m = 5, maxit = 5)
#' 
#' # a single imputation
#' complete1 <- complete(my.mids)
#' 
#' # Backwards elimination for fixed effect models
#' ForwardSelect(y ~ x + w + z, data = complete1)
#' ForwardSelect(y ~ x + w + z, data = my.mids)
#' 
#' # Backwards elimination for mixed (fixed and random) models
#' ForwardSelect(y ~ (1 | factor.1) + x + w + z, data = complete1)
#' ForwardSelect(y ~ (1 | factor.1) + x + w + z, data = my.mids)
#' 
#' @references
#' Douglas Bates and Martin Maechler (2010).
#' lme4: Linear mixed-effects models using S4 classes. R package
#'  version 0.999375-37. http://CRAN.R-project.org/package=lme4
#' 
#' Stef van Buuren, Karin Groothuis-Oudshoorn (2011).
#' mice: Multivariate Imputation by Chained Equations in R.
#' Journal of Statistical Software, 45(3), 1-67.
#' URL http://www.jstatsoft.org/v45/i03/.
#' 
#' @export
ForwardSelect <- function(formula, data, cutoff = .05, family = gaussian,
                          ts.model = NULL, verbose = TRUE) {
  forward.model <- SequentiallyBuildModel(formula, data, cutoff, family,
                                          ts.model, type = "forward", verbose)
  return(forward.model)
}

SequentiallyBuildModel <- function(formula, data, cutoff = .05,
                                   family = gaussian, ts.model, type,
                                   verbose = FALSE) {
  # Private function encaptulating the common logic of the stepwise procedures
  #
  # Args:
  #  formula: a formula with optionally contains random effects
  #  data: a data.frame or mids object
  #  cutoff: the alpha-to-stay or alpha-to-enterparameter
  #  family: the family functions in R
  #  type: either "forward" or "backward"
  #  ts.model:
  #
  # Returns: a gfo model object
  all.terms <- attributes(terms(formula))$term.labels
  fixed.terms <- all.terms[!grepl("\\|", all.terms)]
  random.terms <- all.terms[grepl("\\|", all.terms)]
  response <- all.vars(formula)[1]

  null.formula <- CreateFormula(response, 1)
  complete1.data <- complete(data, 1)
  null.model <- glm(formula = null.formula, family = family,
                    data = complete1.data)

  iteration.terms <- fixed.terms
  term.coef.map <- list()
  history <- list()
  if (type == "forward") {
    weak.term.queue <- c() # For backward elimination
    model.terms <- character(0)
    history[[1]] <- list(formula = CreateFormula(response, 1), fsr.est = 0,
                         term = "")
  }
  if (type == "backward") {
    model.terms <- fixed.terms
    history[[1]] <- list(formula = CreateFormula(response, iteration.terms,
                                                 random.terms),
                         fsr.est = 0, term = "")
  }
  continue <- TRUE
  i <- 0
  while (continue) {
    i <- i + 1
    if (type == "backward") {
      form <- CreateFormula(response, iteration.terms, random.terms)
      test.model <- current.model <-
        GetEstimates(data, form, family, null.model, random.terms, ts.model)
    }
    p.values <- c()
    if (verbose) {
      cat("\n  **** Stepwise Iteration", i, "****\n")
      cat("Testing individual terms...\n--------\n")
    }
    for (var in iteration.terms) {
      form <- CreateFormula(response, union(model.terms, var), random.terms)
      test.model <- GetEstimates(data, form, family, null.model, random.terms,
                                 ts.model)
      test.coefs <- GetCoefNames(var, complete(data))
      new.pvalue <- tryCatch(GetWaldPValue(test.model, test.coefs),
                             error = function(e) {
                               ifelse(type == "backward", -999, 999)})
      term.coef.map[[var]] <- test.coefs
      p.values <- c(p.values, new.pvalue)
      if (verbose) {
        cat("var: ", var, ", pvalue: ", sprintf("%.4f", new.pvalue), " \n")
      }
    }
    iter.list <- RunSelectionCore(type, i, iteration.terms, p.values,
                                           weak.term.queue, cutoff, verbose,
                                           fixed.terms, model.terms, continue)

    form <- CreateFormula(response, iter.list$model.terms, random.terms)
    iteration.terms <- iter.list$iteration.terms
    model.terms <- iter.list$model.terms
    weak.term.queue <- iter.list$weak.term.queue # Backward elimination only
    p.values <- iter.list$p.values
    continue <- iter.list$continue
    if (iter.list$term.added.or.dropped) {
      history[[length(history) + 1]] <- list(formula = form,
                                             fsr = iter.list$fsr,
                                             term = iter.list$term.of.interest)
    }
    if (verbose) {
      cat("\n\nIteration", i, "formula:\n", deparse(form), "\n")
    }
  }
  current.model <- GetEstimates(data, form, family, null.model, random.terms,
                                ts.model)
  current.model$p.values <- p.values
  current.model$term.coef.map <- term.coef.map
  current.model$initial.fixed.terms <- fixed.terms
  current.model$final.fixed.terms <- model.terms
  current.model$random.terms <- random.terms
  current.model$data <- data
  current.model$call.formula <- formula
  current.model$var.select.type <- type
  current.model$var.select.cutoff <- cutoff
  current.model$history <- history
  print(current.model)
  return(current.model)
}

RunSelectionCore <- function(var.select.type, iteration, iteration.terms,
                             p.values,
                             weak.term.queue, cutoff, verbose,
                             fixed.terms, model.terms, continue) {
  if (var.select.type == "forward") {
    return(ForwardSelectCore(iteration, iteration.terms, p.values,
                             weak.term.queue, cutoff, verbose,
                             fixed.terms, model.terms, continue))
  }
  if (var.select.type == "backward") {
    return(BackwardEliminationCore(iteration, iteration.terms, p.values,
                             weak.term.queue, cutoff, verbose,
                             fixed.terms, model.terms, continue))
  }
}

ForwardSelectCore <- function(iteration, iteration.terms, p.values,
                              weak.term.queue, cutoff, verbose,
                              fixed.terms, model.terms, continue) {
  # The stepwise logic unique to forward selection
  msv.index <- which.min(p.values)
  smallest.p.value <- p.values[msv.index]
  most.sig.var <- iteration.terms[msv.index]
  fsr.est <- 0 # if you never select anything, you have no false discoveries
  term.added.or.dropped <- FALSE
  if (length(iteration.terms) > 0 && smallest.p.value <= cutoff) {
    model.terms <- union(model.terms, most.sig.var)
    iteration.terms <- iteration.terms[-msv.index]
    term.added.or.dropped <- TRUE
    if (verbose) {
      cat("Forward selection iteration", iteration, ",adding term",
          most.sig.var, "to model\n")
    }
    fsr.est <- GetFastFSR(n.total.vars = length(fixed.terms),
                          n.model.vars = length(model.terms),
                          alpha = smallest.p.value)
  }
  if (length(iteration.terms) == 0 || smallest.p.value > cutoff) {
     continue <- FALSE
     if (verbose) cat("No further terms have p-values <", cutoff, "\n")
  }
  return(list(fsr = fsr.est, model.terms = model.terms,
              iteration.terms = iteration.terms,
              weak.term.queue = weak.term.queue, p.values = p.values,
              continue = continue, term.of.interest = most.sig.var,
              term.added.or.dropped = term.added.or.dropped))
}

BackwardEliminationCore <- function(iteration, iteration.terms, p.values,
                                    weak.term.queue, cutoff, verbose,
                                    fixed.terms, model.terms, continue) {
  # The stepwise logic unique to backward selection
    term.added.or.dropped <- FALSE
    lsv.index <- which.max(p.values)
    if (max(p.values) == -999) {  # Bad run. Drop term from weak.term.queue.
      least.sig.var <- weak.term.queue[1]
      weak.term.queue <- weak.term.queue[-1]
      largest.p.value <- 1
      warning("All runs failed. Using p-values from previous runs")
    } else {
      largest.p.value <- p.values[lsv.index]
      # The smallest alpha such that all remaining  terms stay in the model
      next.largest.p.value <- ifelse(length(p.values) > 1,
                                     max(p.values[-lsv.index]),
                                     largest.p.value)
      least.sig.var <- iteration.terms[lsv.index]
      # creation of weak.term.queue: keep a sequence of the least significant
      # variables from a successful run. Use it when an iteration fails.
      weak.term.queue <- iteration.terms[order(p.values, decreasing = TRUE)][-1]
    }
    if (length(iteration.terms) > 0 && largest.p.value > cutoff) {
      p.values <- p.values[-lsv.index]
      names(p.values) <- iteration.terms[-lsv.index]
      model.terms <- iteration.terms <- setdiff(iteration.terms, least.sig.var)
      if (verbose) {
        term.added.or.dropped <- TRUE
        cat("Backward Elimination iteration", iteration, ", dropping term",
            least.sig.var, "from model\n")
      }
    }
    if (length(iteration.terms) == 0 || largest.p.value <= cutoff) {
      continue <- FALSE
      model.terms <- iteration.terms
      names(p.values) <- iteration.terms
      if (verbose) cat("All remaining terms have p-values <", cutoff, "\n")
    }
    fsr.est <- GetFastFSR(n.total.vars = length(fixed.terms),
                          n.model.vars = length(model.terms),
                          alpha = next.largest.p.value)

    return(list(fsr = fsr.est, model.terms = model.terms,
                iteration.terms = iteration.terms,
                weak.term.queue = weak.term.queue, p.values = p.values,
                continue = continue, term.of.interest = least.sig.var,
                term.added.or.dropped = term.added.or.dropped))
}
