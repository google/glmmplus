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

# A description of the following data generating process can be found here:
# http://goo.gl/eIIxg0

# TODO (baogorek): Add test for AR1 parameter

N <- 500
n <- 15
ar <- 0.8
sd.epsilon <- 0.5

dat.list <- list()
ids.list <- list()
fixed.list <- list()
ran.effs <- .75 * rnorm(N)
x <- rnorm(N)
for (i in 1:N) {
  dat.list[[i]] <- 0 + ran.effs[i] + .7 * x[i] +
                   arima.sim(n = n, list(ar = ar, sd = sd.epsilon))
  ids.list[[i]] <- rep(i, n)
  fixed.list[[i]] <- rep(x[i], n)
}

longit.df <- data.frame(id = unlist(ids.list),
                        x = unlist(fixed.list),
                        y = unlist(dat.list),
                        y.binary = as.numeric(unlist(dat.list) < 0),
                        subject = sample(1:15, N*n, replace = TRUE)
)

longit.df[sample(1:N, round(N / 4)), "x"] <- NA
library(mice)
mids <- mice(longit.df)

library(nlme)
my.nlme <- lme(y ~ x, data = complete(mids), random = ~ 1 | id,
               cor = corAR1(0.5, form = ~ 1 | id))
library(lme4)
my.lmer <- lmer(y ~ (1 | id) + x, data = complete(mids))

TestNLMEmatch <- function() {
  my.gfo <- FitModel(y ~ (1 | id) + x, data = complete(mids),
                     ts.model = "ar1")
  checkEquals(round(coef(my.gfo)[2], 3), round(fixef(my.nlme)[2], 3))
}

TestNonGaussianException <- function() {
  checkException(FitModel(y.binary ~ (1 | id) + x, data = complete(mids),
                          family = binomial, ts.model = "ar1"))

  checkException(FitModel(y.binary ~ (1 | id) + x, data = mids,
                          family = binomial, ts.model = "ar1"))

  # Look for "Sorry, time series models have not been implemented
  # for non-gaussian models"
}

TestMidsTsModel <- function() {
  gfo.mids <- FitModel(y ~ (1 | id) + x, data = mids,
                     ts.model = "ar1")
  checkEquals(round(coef(gfo.mids)[2], 1), round(fixef(my.nlme)[2], 1))
}

