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

WideToLong <- function(object, ...) UseMethod('WideToLong', object)

# Removes leading and trailing white space
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

tr <- function(mat) {
  # Creates a matrix trace
  #
  # Args:
  #  mat: a matrix
  #
  # Returns: a scalar trace
  trace <- sum(diag(1, nrow(mat), nrow(mat)) * mat)
  return(trace)
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
      test.coefs <- paste0(var, colnames(contrasts(obj)))
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

# Not sure if this is necessary
RemoveRandomEffectVars <- function(var.vec) {
 # Removes random effect terms from a term vector
 #
 # Args:
 #  var.vec: a character vector of random terms
 #
 # Returns: a vector of variable names
 non.random.terms <- var.vec[!grepl("\\|", var.vec)]
 return(trim(non.random.terms))
}

MakeSplineStr <- function(var.name, data, interior.knots) {
  # Creates an ns function string designed for reuse on future data sets
  #
  # Args:
  #  var.name: the variable name of the independent variable
  #  data: a data.frame
  #  interior.knots: a numeric vector of parametric spline knots
  library(splines)
  str.1 <- paste0("ns(", var.name, ", knots = c(")
  str.2 <- paste(interior.knots, collapse = ", ")
  str.3 <- "), Boundary.knots = c("
  str.4 <- paste(round(range(data[, var.name]), .2), collapse = ", ")
  return(paste0(str.1, str.2, str.3, str.4, "))"))
}

# TODO: Make this work for more complicated random effect situations
CreateFormula <- function(response, predictors, random.terms = character(0)) {
  # Creates a formula based on terms passed in as arguments
  #
  # Args:
  #  reponse: the response variable's name
  #  predictors: a character vector of predictor names
  #  random.terms: a character vector of random terms as they are presented in
  #                terms(formula)
  #
  # Returns: a formula
  if (length(predictors) == 0) predictors <- "1"
  random.part <- ifelse(length(random.terms) == 0, "",
                        paste0("(", paste(random.terms, collapse = ") + ("),
                               ") + "))
  fixed.part <- paste(predictors, collapse = "+")
  formula.str <- paste(response, "~", random.part, fixed.part)
  return(formula(formula.str))
}

WideToLong.data.frame <- function(data, id.name, response.base,
                                  time.varying.bases, sep = ".") {

  response.regex <- paste0("^", gsub("\\.", "\\\\.",
                           paste0(response.base, sep)), "(.*)$")
  response.names <- grep(response.regex, names(data), value = TRUE)
  period.names <- gsub(response.regex, "\\1", response.names)

  
  wide <- data[, setdiff(names(data), c(response.names))] 
  # putting period first guarantees it will keep its order
  long <- expand.grid(period.names, data[, id.name]) [, c(2, 1)]
  names(long) <- c(id.name, "period")

  for (var.base in c(response.base, time.varying.bases)) {
    # Get the names of the variables
    # knock them out of wide
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

  #static.expanded <- merge(wide, data.frame(period = period.names))

  long.df <- merge(long, wide)
  return (long.df)
}

