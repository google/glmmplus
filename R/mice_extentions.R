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

ImputeData <- function(data, m = 10, maxit = 15, droplist = NULL) {
  if (length(intersect(names(data), droplist)) < length(droplist)) {
    stop("Droplist variables not found in data set")
  }
  require(mice)
  predictorMatrix = (1 - diag(1, ncol(data)))
  for (term in droplist) {
    drop.index <- which(names(data) == term)
    predictorMatrix[, drop.index] <- 0
  }
  mids.out <- mice(data, m = m, maxit = maxit,
                   predictorMatrix = predictorMatrix)
  return(mids.out)
}

WideToLong <- function(object, ...) UseMethod('WideToLong', object)
WideToLong.mids <- function(wide.mids, id.name, response.base,
                            time.varying.bases = NULL, sep = ".") {

  new.data <- WideToLong(wide.mids$data, id.name, response.base,
                         time.varying.bases, sep)
  new.imp <- list()
  m <- wide.mids$m  # Number of imputed data sets
  for (var in names(new.data)) {
    na.indices <- which(is.na(new.data[, var]))
    n.mis <- length(na.indices)
    M <- matrix(rep(0, n.mis * m), ncol = m)

    rownames(M) <- na.indices 
    colnames(M) <- 1:m
    # iterate throu the complete sets and fill them in
    for (i in 1:m) {
      complete.i <- WideToLong(complete(wide.mids, i), id.name, response.base,
                               time.varying.bases, sep)

      M[, i] <- complete.i[na.indices, var]
    }
    new.imp[[var]] <- as.data.frame(M)
  }
  long.mids <- list(m = m, imp = new.imp, data = new.data)
  class(long.mids) <- "mids"
  return(long.mids)
}

`[.mids` <- function(obj, first, second) {
  # Row level subsetting of a mids object
  #
  # Args:
  #  obj: a mids object
  #  first: the row subset
  #  second: the column subset (TODO: Not currently implemented)
  #
  # Returns: a mids object
  new.data <- obj$data[first, ]
  new.imp <- list()
  # imp is A list of data frames, one for each variable with missing values
  # Each data frame as m columns, one for each imputations.
  for (var in names(obj$imp)) {
    rows.found <- rownames(obj$imp[[var]]) %in% rownames(new.data)
    if (sum(rows.found) == 0 || is.null(obj$imp[[var]])) {
      new.imp[var] <- list(NULL)
    } else {
      new.imp[[var]] <- obj$imp[[var]][rows.found, ]
    }
  }
  new.obj <- list(m = obj$m, imp = new.imp, data = new.data)
  class(new.obj) <- "mids"
  return(new.obj);
}

# TODO: bring this back. Add in plot(mice).
# Also, if I can't serve it, would rather have pdf and html
#CreateImputationReport <- function(x.df, imputation.fn, out.path){
#  # Create a report of the data set variables and the imputation
#  library(knitr)
#
#  x.df.imp <- imputation.fn(x.df)
#  x.df.complete  <- complete(x.df.imp, 1)
#
#  knit2html("reports/imputations.Rmd")
#  # move the file to the user specified directory
#  file.rename("imputations.html",out.path)
#  # clean up
#  unlink("imputations.md")
#  unlink("figure", recursive=TRUE)
#}

