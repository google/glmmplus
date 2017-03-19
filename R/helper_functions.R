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

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

GetFnArg <- function(fn.str) {
  # Grabs the first argument from a function
  #
  # Args:
  #  fn.str: a string representation of a function
  #
  # Returns: a numeric string representation of the first argument
  s <- gsub("[a-zA-Z._\\s]+\\(\\s?([a-zA-Z._]+)\\s?[\\)|,].*", "\\1", fn.str)
  return(trim(s))
}

RelabelSplits <- function(vec, values.to.shuffle) {
  # Reshuffle a subset of an integer vector
  #
  # Args:
  #  vec: an integer vector
  #  values.to.shuffle: the values in vec that should be reshuffled
  #
  # Returns: a vector with the same values as vec but in a different order
  positions <- vec %in% values.to.shuffle
  values <- vec[positions]
  shuffled.values <- sample(values, length(values))
  new.vec <- vec
  new.vec[positions] <- shuffled.values
  return(new.vec)
}

tr <- function(mat) { # create a matrix trace
  sum(diag(1, nrow(mat), nrow(mat)) * mat)
}

GetCoefNames <- function(var, complete.df) {
  # Gets the coefficient names from a variable
  #
  # Args:
  #  var: a variable name (character)
  #  complete.df: a data.frame
  #
  # Return: a character vector of variable names
  obj <- with(complete.df, eval(parse(text = var)))
  if (grepl("\\(", var) && !is.null(dim(obj))) {  # e.g., an ns() spline fn
    test.coefs <- paste0(var, 1:ncol(obj))
  } else {
      test.coefs <- var
    if (is.factor(obj)) {
      test.coefs <- paste0(var, colnames(contrasts(factor(obj))))
    }
    else if (is.character(obj)) {
      test.coefs <- paste0(var, levels(factor(obj))[-1])
    }
    else if (is.matrix(obj)) {
      test.coefs <- paste0(var, colnames(obj))
    }
  }
  return(test.coefs)
}

ReformatNames <- function(old.names) {
  # Computes variable names that in accordance with the Google style guide
  #
  # Args:
  #   old.names: a vector of variable names
  #
  # Returns:
  #   A vector of names that have been transformed
  new.names <- tolower(old.names)
  new.names <- gsub("_", ".", new.names)
  return(new.names)
}

IsContinuous <- function(data, min.levels = 5) {
  # Finds continuous variables (in a practical sense) from a data frame
  #
  # Args:
  #   data: a data frame
  #   min.levels: The minimum number of unique values that a variable must
  #               have to be considered "continuous"
  #
  # Returns:
  #   A logical vector where TRUE corresponds to continuous
  is.cont <- sapply(1:ncol(data), FUN = function(i) {is.numeric(data[, i]) &&
        length(unique(data[, i])) > min.levels})
  return(is.cont)
}

GetStandardDevs <- function(data) {
  # Get standard deviations of all continuous variables in a data set
  #
  # Args:
  #   data: a data.frame
  #
  # Returns: A numeric vector of standard deviations
  sd.vec <- sapply(data[, IsContinuous(data)], function(x) sd(x, na.rm = TRUE))
  return(sd.vec)
}

GetRelatedVars <- function(base, ref.df) {
  # Finds all variables with a given base
  #
  # Args:
  #   base: a character prefix
  #   ref.df: a reference data frame
  #
  # Returns: a character vector of variable names
    pattern <- paste0("^", base)
    occurances <- grepl(pattern, names(ref.df))
    variables <- names(ref.df)[occurances]
    return(variables)
}

strip.rnd.code <- function(rnd.str) {
  # Extracts the random term from an lmer parentheses block
  #
  # Args:
  #  rnd.str: a piece of a mer formula containing a random term
  #
  # Returns: a variable name
  return(trim(sub("\\([^\\|]+\\|\\s*(.*)\\s*\\)", "\\1", rnd.str)))
}

#TODO: Make work for nested and more complicated cases
GetRandomEffectVars <- function(var.vec) {
  # Extracts the variables from a vector of random formula terms
  #
  # Args:
  #  var.vec: a character vector of random terms
  #
  # Returns: a vector of variable names
   random.terms <- var.vec[grepl("\\|", var.vec)]
   return.vec=c()
   for (i in 1:length(random.terms)) {
   return.vec <- c(return.vec, strip.rnd.code(random.terms[i]))
   }
   return(trim(return.vec))
}

RemoveRandomEffectVars <- function(var.vec) {
 # Removes random effect terms from a term vector
 non.random.terms <- var.vec[!grepl("\\|", var.vec)]
 return(trim(non.random.terms))
}

#' Make a spline formula string corresponding to unambiguous spline
#'
#' The ns() function is convenient, but without specifiying the knots and the
#' boundary knots, the information neccessary to reproduce the function is in
#' the data set itself. This function creates a string that can be used in
#' future formulas and does not depend on the original data set for
#' reproducibility.
#'
#' @param var.name name of variable for use with spline equation
#' @param data data set that will be used to construct the spline string
#' @param interior.knots the interior knots of the spline
#'
#' @export
MakeSplineStr <- function(var.name, data, interior.knots) {
  str.1 <- paste0("ns(", var.name, ", knots = c(")
  str.2 <- paste(interior.knots, collapse = ", ")
  str.3 <- "), Boundary.knots = c("
  str.4 <- paste(round(range(data[, var.name]), .2), collapse = ", ")
  return(paste0(str.1, str.2, str.3, str.4, "))"))
}

#' Formula creation from string arguments
#'
#' @param response the response variable's name
#' @param predictors a character vector of predictor names
#' @param random.terms a character vector of random terms as they are presented in
#'                terms(formula)
#'
#' @export
CreateFormula <- function(response, predictors, random.terms = character(0)) {
  if (length(predictors) == 0) {
    predictors <- "1"
  }
  random.part <- ifelse(length(random.terms) == 0, "",
                        paste0("(", paste(random.terms, collapse = ") + ("),
                               ") + "))
  fixed.part <- paste(predictors, collapse = "+")
  formula.str <- paste(response, "~", random.part, fixed.part)
  new.formula <- as.formula(formula.str)
  return(new.formula)
}

#' WideToLong: Converting from wide to long formats
#'
#' In longitudinal or other multiple response studies, data presented in a long
#' format will often feature dependence between rows. While this is the
#' preferred format for lme4, such a format would hide important information
#' from multiple imputation models and make the MAR assumption less plausible.
#' Hense, the suggestion is to impute data in a wide format, where rows are
#' again independent, and then return the mids object to a long format for use
#' with FitModel, ForwardSelect, or BackwardEliminate.
#'
#'
#' @param data A data frame or mids object in "wide" format. Specifically, both
#'             the response and any time varying covariates should be specified
#'             as multiple columns with the same base name, but a different
#'             suffix. The suffix values will be the future period labels.
#'  
#' @param id.name The name of the identifying variable, a character string.
#'  
#' @param response.base The common prefix for the response variable, a character
#'                 string.
#' @param time.varying.bases A character vector of name prefixes for
#'        time-varying covariates.
#' @param sep The character delimiter separating the variable name base from
#'        the period identifier.
#' 
#' @seealso \code{\link{LongToWide}}
#' @examples
#' wide.df <- data.frame(pid           = 1:100,
#'                       my.response.1 = rnorm(100),
#'                       my.response.2 = rnorm(100),
#'                       x.1           = rnorm(100),
#'                       x.2           = rnorm(100))
#' # add missingness
#' 
#' wide.df[25:50, "my.response.2"] <- NA
#' wide.df[45:55, "x.1"] <- NA
#' 
#' wide.mids <- ImputeData(wide.df, droplist = c("pid"))
#' long.mids <- WideToLong(wide.mids, "pid", "my.response", c("x"), sep = ".")
#' 
#' my.model <- FitModel(my.response ~ (1 | pid) + x, data = long.mids)
#' summary(my.model)
#' @references
#' Stef van Buuren, Karin Groothuis-Oudshoorn (2011).
#' mice: Multivariate Imputation by Chained Equations in R. Journal of
#'       Statistical Software, 45(3), 1-67. URL
#'       http://www.jstatsoft.org/v45/i03/.
#'
#' @export
WideToLong <- function(data, id.name, response.base, time.varying.bases, sep) {
  UseMethod('WideToLong', data)
}

#' @rdname WideToLong 
#' @export
WideToLong.data.frame <- function(data, id.name, response.base,
                                  time.varying.bases = NULL, sep = ".") {

  response.regex <- paste0("^", gsub("\\.", "\\\\.",
                           paste0(response.base, sep)), "(.*)$")
  response.names <- grep(response.regex, names(data), value = TRUE)
  period.names <- gsub(response.regex, "\\1", response.names)

  
  wide <- data[, setdiff(names(data), c(response.names))] 
  # putting period first guarantees it will keep its order
  long <- expand.grid(period.names, data[, id.name]) [, c(2, 1)]
  names(long) <- c(id.name, "period")

  for (var.base in c(response.base, time.varying.bases)) {
    # Get the names of the variables and take them out of wide
    regex <- paste0("^", gsub("\\.", "\\\\.",
                    paste0(var.base, sep)), "(.*)$")
    names <- grep(regex, names(data), value = TRUE)
    subset <- data[, names] 
    ts.vec <- unlist(as.data.frame(t(subset)))
    long[, var.base] <- ts.vec
    for (name in names) {
      wide[, name] <- NULL
    }
  }
  long.df <- merge(long, wide)
  return (long.df)
}

#' LongToWide: Convert nested long structures to wide multivariate structures
#'
#' In longitudinal or other multiple response studies, data presented in a long
#' format will often feature dependence between rows. While this is the
#' preferred format for lme4, such a format would hide important information
#' from multiple imputation models and make the MAR assumption less plausible.
#' Hense, the suggestion is to impute data in a wide format, where rows are
#' again independent and then return the mids object to a long format for use
#' with FitModel, ForwardSelect, or BackwardEliminate.
#' 
#' @param data A data frame or mids object in "long" format owing to multiple
#'             measurements within the same subject.
#' @param id.name The subject id, a character string.
#' @param period.name The repeated measurement (within subject) identifier.
#'                    In a longitudinal study, this will be time.
#' @param time.varying.vars A character vector of variable names that take
#'        multiple values per subject (in different rows)
#' @param  sep The character delimiter by which to separate the variable name
#'             base from the period identifier.
#' 
#' @seealso \code{\link{WideToLong}}
#' 
#' @examples
#' # Example of the long-to-wide, impute, wide-to-long strategy
#' library(glmmplus)
#' data(nls.97)
#' nls.97[1:10, 1:4]
#' 
#' nls.wide <- LongToWide(nls.97, id.name = "PUBID.1997", period.name = "age",
#'                        time.varying.vars = c("math.cs"))
#' 
#' nls.wide[1:2, c(1:2, 20:29)]
#' mids <- ImputeData(nls.wide, m = 5, maxit = 15, droplist = c("PUBID.1997"))
#' mids.long <- WideToLong(mids, "PUBID.1997", "math.cs")
#' 
#' @references
#' Stef van Buuren, Karin Groothuis-Oudshoorn (2011).
#' mice: Multivariate Imputation by Chained Equations in R.
#'       Journal of Statistical Software, 45(3), 1-67. URL
#'       http://www.jstatsoft.org/v45/i03/.
#' @export
LongToWide <- function(data, id.name, period.name, time.varying.vars, sep) {
  UseMethod('LongToWide', data)
}

#' @rdname LongToWide
#' @export
LongToWide.data.frame <- function(data, id.name, period.name,
                                  time.varying.vars, sep = ".") {
  data <- data[order(data[, id.name], data[, period.name]), ]
  stable.vars <- setdiff(names(data), c(period.name, time.varying.vars))
  if (length(stable.vars) == 1) {
    time.stable.df <- data.frame(unique(data[, id.name]))
    names(time.stable.df) <- id.name
  } else {
    time.stable.df <- unique(data[, stable.vars])
  }
  period.names <- unique(data[, period.name])
  wide.df <- time.stable.df
  for (var in time.varying.vars) {
    for (period in period.names) {
     period.df <- data[data[, period.name] == period, c(id.name, var)]
     names(period.df) <-  c(id.name, paste(var, period, sep = sep))
     wide.df <- merge(wide.df, period.df, by = id.name, all.x = TRUE)
    }
  }
  return(wide.df)
}
