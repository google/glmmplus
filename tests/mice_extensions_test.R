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

wide.df <- data.frame(pid           = 1:100,
                      my.response.1 = rnorm(100),
                      my.response.2 = rnorm(100),
                      x.1           = rnorm(100),
                      x.2           = rnorm(100))

wide.df[25:50, "my.response.2"] <- NA
wide.df[45:55, "x.1"] <- NA

sink("/dev/null")
mids <- ImputeData(wide.df, droplist = c("pid"))
sink()

TestImputationWrapper <- function() {
   
  checkException(ImputeData(wide.df, droplist = c("not_here")))

  imp.predictors <- apply(mids$predictorMatrix, 2, FUN = sum)
  used.predictors <- sort(names(imp.predictors[imp.predictors > 0]))
  checkEquals(used.predictors, setdiff(names(wide.df), "pid"))
}

TestMidsSubset <- function() {
  subset.mids <- mids[50:100, ]
  checkTrue(all(complete(mids, 1)[50:100, ] == complete(subset.mids, 1)))
  checkTrue(all(complete(mids, 2)[50:100, ] == complete(subset.mids, 2)))
}

TestWideToLong <- function() {
  long.mids <- WideToLong(mids, "pid", "my.response", c("x"), sep = ".")
  long.df <- complete(long.mids, 1)
  
  checkEquals(nrow(long.df), 2 * nrow(wide.df))
  checkEquals(sort(names(long.df)),
              sort(c("pid", "my.response", "period", "x")))

  # checking cases with no missing values
  checkEquals(wide.df[1, "my.response.1"], long.df[1, "my.response"])
  checkEquals(as.character(long.df[1, "period"]), "1")

  checkEquals(wide.df[1, "my.response.2"], long.df[2, "my.response"])
  checkEquals(as.character(long.df[2, "period"]), "2")

  checkEquals(wide.df[1, "x.1"], long.df[1, "x"])
  checkEquals(wide.df[1, "x.2"], long.df[2, "x"])

  checkTrue(!all(complete(long.mids, 1) == complete(long.mids, 2)))
}
