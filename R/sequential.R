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

BackwardEliminate <- function(formula, data, cutoff = .05,
                              family = gaussian, ts.model = NULL,
                              verbose = TRUE) {
  # Model selection by backward elimination from a full model to a reduced one
  #
  # Args:
  #  formula: a formula with optionally contains random effects
  #  data: a data.frame or mids object
  #  cutoff: the alpha-to-stay parameter
  #  family: the family functions in R
  #  verbose: logical - print out results from each iteration?
  #
  # Returns: a gfo model object
  final.model <- SequentiallyBuildModel(formula, data, cutoff, family,
                                        ts.model, type = "backward", verbose)
  return(final.model)
}

ForwardSelect <- function(formula, data, cutoff = .05, family = gaussian,
                          ts.model = NULL, verbose = TRUE) {
  # Model selection by forward selection from a null model onwards
  #
  # Args:
  #  formula: a formula with optionally contains random effects
  #  data: a data.frame or mids object
  #  cutoff: the alpha-to-enter parameter
  #  family: the family functions in R
  #
  # Returns: a gfo model object
  forward.model <- SequentiallyBuildModel(formula, data, cutoff, family,
                                          ts.model, type = "forward", verbose)
  if (length(forward.model$final.fixed.terms) == 0) {
    final.model <- forward.model
  } else {
    final.model <- FitModel(forward.model$formula, data, family)
    # Replacing data that was rewritten above
    final.model$var.select.type <- "forward"
    final.model$var.select.cutoff <- forward.model$var.select.cutoff
    final.model$fsr.est <- forward.model$fsr.est
  }
  return(final.model)
}

SequentiallyBuildModel <- function(formula, data, cutoff = .05,
                                   family = gaussian, ts.model, type, verbose) {
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
  library(splines)
  all.terms <- attributes(terms(formula))$term.labels
  fixed.terms <- all.terms[!grepl("\\|", all.terms)]
  random.terms <- all.terms[grepl("\\|", all.terms)]
  response <- all.vars(formula)[1]
  null.model <- glm(CreateFormula(response, 1), family, complete(data))
  iteration.terms <- fixed.terms
  term.coef.map <- list()
  if (type == "forward") model.terms <- character(0)
  if (type == "backward") model.terms <- fixed.terms
  continue <- TRUE
  i <- 0
  while (continue) {
    i <- i + 1
    if (type == "backward") {
      form <- CreateFormula(response, iteration.terms, random.terms)
      test.model <- current.model <- GetEstimates(data, form, family,
                                                  null.model, random.terms,
                                                  ts.model)
    }
    p.values <- c()
    if (verbose) {
    cat("\n  **** Stepwise Iteration", i, "****\n")
    cat("Testing individual terms...\n--------\n")
    }
    for (var in iteration.terms) {
      if (type == "forward") {  # fit model for each variable
        form <- CreateFormula(response, union(model.terms, var), random.terms)
        test.model <- GetEstimates(data, form, family, null.model, random.terms,
                                   ts.model)
      }
      test.coefs <- GetCoefNames(var, complete(data))
      new.pvalue <- tryCatch(GetWaldPValue(test.model, test.coefs),
                             error = function(e) {
                               bad.p <- 999
                               if (type == "backward") bad.p <- -999
                               return(bad.p)
                             })
      term.coef.map[[var]] <- test.coefs
      if (new.pvalue < min(c(1, p.values)) && type == "forward") {
        current.model <- test.model
      }
      p.values <- c(p.values, new.pvalue)
      if (verbose) {
        cat("var: ", var, ", pvalue: ", sprintf("%.4f", new.pvalue), " \n")
      }
    }
    if (verbose) cat("--------\n")
    if (type == "forward") ForwardSelectCore()
    if (type == "backward") BackwardEliminationCore()
    form <- CreateFormula(response, model.terms, random.terms)
    if (verbose) cat("Iteration", i, "formula:\n", deparse(form), "\n")
  }
  current.model <- GetEstimates(data, form, family, null.model, random.terms,
                                ts.model)
  names(p.values) <- iteration.terms

  fsr.est <- GetFastFSR(n.total.vars = length(fixed.terms),
                        n.model.vars = length(model.terms),
                        alpha = cutoff)
  if (verbose) cat("Final FSR estimate : ", round(fsr.est, 3), "\n")

  current.model$p.values <- p.values
  current.model$term.coef.map <- term.coef.map
  current.model$initial.fixed.terms <- fixed.terms
  current.model$final.fixed.terms <- model.terms
  current.model$random.terms <- random.terms
  current.model$data <- data
  current.model$call.formula <- formula
  current.model$var.select.type <- type
  current.model$var.select.cutoff <- cutoff
  current.model$fsr.est <- fsr.est
  return(current.model)
}

ForwardSelectCore <- function() {
  # The stepwise logic unique to forward selection
  with(parent.frame(), {
  msv.index <- which.min(p.values)
  smallest.p.value <- p.values[msv.index]
  most.sig.var <- iteration.terms[msv.index]

  if (length(iteration.terms) > 0 && smallest.p.value <= cutoff) {
    model.terms <- union(model.terms, most.sig.var)
    iteration.terms <- iteration.terms[-msv.index]
    if (verbose) {
      cat("Forward selection iteration", i, ",adding term", most.sig.var,
          "to model\n")
    }
    fsr.est <- GetFastFSR(n.total.vars = length(fixed.terms),
                          n.model.vars = length(model.terms),
                          alpha = smallest.p.value)
    if (verbose) cat("\n** FSR estimate : ", round(fsr.est, 3), "***\n\n")

  }
  if (length(iteration.terms) == 0 || smallest.p.value > cutoff) {
     continue <- FALSE
     if (verbose) cat("No further terms have p-values <", cutoff, "\n")
  }
  })
}

BackwardEliminationCore <- function() {
  # The stepwise logic unique to backward selection
  with(parent.frame(), {
    lsv.index <- which.max(p.values)
    if (max(p.values) == -999) {  # Bad run. Drop term from weak.term.queue.
      least.sig.var <- weak.term.queue[1]
      weak.term.queue <- weak.term.queue[-1]
      largest.p.value <- 1
      warning("All runs failed. Using p-values from previous runs")
    } else {
      largest.p.value <- p.values[lsv.index]
      least.sig.var <- iteration.terms[lsv.index]
      # creation of weak.term.queue: keep a sequence of the least significant
      # variables from a successful run. Use it when an iteration fails.
      weak.term.queue <- iteration.terms[order(p.values, decreasing = TRUE)][-1]
    }
    if (length(iteration.terms) > 0 && largest.p.value > cutoff) {
      model.terms <- iteration.terms <- setdiff(iteration.terms, least.sig.var)
      if (verbose) {
        cat("Backward Elimination iteration", i, ", dropping term", least.sig.var,
            "from model\n")
      }

    second.largest.p.value = max(p.values[-largest.p.value])
    fsr.est <- GetFastFSR(n.total.vars = length(fixed.terms),
                          n.model.vars = length(model.terms),
                          alpha = second.largest.p.value)
    if (verbose) cat("\n***FSR estimate : ", round(fsr.est, 3), "****\n\n")
    }

    if (length(iteration.terms) == 0 || largest.p.value <= cutoff) {
      continue <- FALSE
      if (verbose) cat("All remaining terms have p-values <", cutoff, "\n")

    }
  })
}
