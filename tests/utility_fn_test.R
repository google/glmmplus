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

TestTrim <- function() {
  test.str1 <- " abc "
  checkEquals(glmmplus:::trim(test.str1), "abc")
}

TestGetFnArg <- function() {
  test.str1 <- "ns( my.var , other.arg = 2)"
  test.str2 <- "ns(my.var , other.arg = 2)"
  test.str3 <- "ns( my.var, other.arg = 2)"
  test.str4 <- "ns(my.var, other.arg = 2)"
  test.str5 <- "cut(x, c(-10, -1, 0, 1, 10))"

  checkTrue(glmmplus:::GetFnArg(test.str1) == "my.var")
  checkTrue(glmmplus:::GetFnArg(test.str2) == "my.var")
  checkTrue(glmmplus:::GetFnArg(test.str3) == "my.var")
  checkTrue(glmmplus:::GetFnArg(test.str4) == "my.var")
  checkTrue(glmmplus:::GetFnArg(test.str5) == "x")
}

TestRelabelSplits <- function() {
  original <- sample(1:5, 1000, replace = TRUE)
  new <- glmmplus:::RelabelSplits(original, 4:5)

  checkTrue(all(new[new %in% 1:3] == original[original %in% 1:3]))
  checkTrue(!all(new[new %in% 4:5] == original[original %in% 4:5]))
  checkTrue(all(table(new) == table(original)))
}

TestContinuousVarOps <- function() {
  test.df <- data.frame(a = integer(7), b = double(7), c = factor(1:7),
                        d = character(7), e = integer(7), f = logical(7))
  test.df$a <- 1:7
  test.df$b <- rnorm(7)
  test.df$d <- as.character(1:7)
  test.df$e <- c(1, 2, 3, 1, 2, 3, 3)
  test.df$f <- c(T, F, T, T, F, T, T)

  logical.vec <- glmmplus:::IsContinuous(test.df)
  checkEquals(names(test.df[, logical.vec]), c("a", "b"))

  std.devs <- glmmplus:::GetStandardDevs(test.df)
  checkEquals(names(std.devs), c("a", "b"))
  checkEquals(as.numeric(std.devs["a"]), sd(1:7))
}

TestMatrixTrace <- function() {
  # Used in GetWaldPValues()
  matrix.1 <- matrix(c(1, 4, 2, 3), ncol = 2)
  matrix.2 <- matrix(c(3, 2, 5, 0, 0, 1, 2, 4, 5), ncol = 3)
  checkEquals(glmmplus:::tr(matrix.1), 4)
  checkEquals(glmmplus:::tr(matrix.2), 8)
}

TestExtractCoefNames <- function() {
  GetCoefNames <- glmmplus:::GetCoefNames
  y <- rnorm(100)
  x <- rnorm(100)
  z <- sample(c("red", "green", "blue"), 100, replace=TRUE)
  w <- cut(x, c(-Inf, -2, 0, 2, Inf))
  my.df <- data.frame(x = x, z = z, w = w)
  my.spline <- MakeSplineStr("x", my.df, c(-1, 1))
  checkEquals(GetCoefNames("x", my.df), "x")
  checkEquals(GetCoefNames("z", my.df), c("zgreen", "zred"))
  checkEquals(GetCoefNames("w", my.df), c("w(-2,0]", "w(0,2]", "w(2, Inf]"))
  checkEquals(GetCoefNames(my.spline, my.df), paste0(my.spline, c(1:3)))
}

TestExtractRandomTerms <- function() {
  random.terms <- glmmplus:::GetRandomEffectVars(c(" (1 | x)", " ( 1 | y)"))
  checkTrue(all(random.terms == c("x", "y")))
}

TestRemoveRandomEffectVars <- function() {
  effects.vec <- c(" (1 | x)", " ( 1 | y)", " z ")
  fixed.terms <- glmmplus:::RemoveRandomEffectVars(effects.vec)
  checkTrue(all(fixed.terms == c("z")))
}

TestReformatNames <- function() {
  test.df <- data.frame(Meas_1 = c(1, 1, 1), Meas_2 = c(2, 1, 2))
  names(test.df) <- ReformatNames(names(test.df))
  checkEquals(names(test.df), c("meas.1", "meas.2"))
}

TestGetRelatedVars <- function () {
  test.df <- data.frame(art.painting = c(1, 1), art.drawing = c(1, 2),
                        prevyr.1 = c(3, 2), prevyr.2 = c(2, 3))

  art.vec <- GetRelatedVars("art.", test.df)
  checkEquals(art.vec, c("art.painting", "art.drawing"))

  prevyr.vec <- GetRelatedVars("prevyr.", test.df)
  checkEquals(prevyr.vec, c("prevyr.1", "prevyr.2"))
}

TestMakeSplineString <- function() {
  my.df <- data.frame(x = c(1, 2, 3))
  my.ns <- MakeSplineStr("x", my.df, c(1.5, 2.5))
  checkEquals(my.ns, "ns(x, knots = c(1.5, 2.5), Boundary.knots = c(1, 3))")
}

TestCreateFormula <- function() {
  form <- glmmplus:::CreateFormula("y", c("x", "z"), c(" 1 | w "))
  checkEquals(form, formula("y ~ (1 | w) + x + z"))
}

TestWideToLongDf <- function() {
  y <- rnorm(4)
  x.t <- rnorm(4) # time varying covariate
  z <- c(3, 2)
  wide.df <- data.frame(p.id = c(1, 2), y.1 = y[c(1, 3)], y.2 = y[c(2, 4)],
                        x.t.1 = x.t[c(1, 3)], x.t.2 = x.t[c(2, 4)], z = z)
  
  long.df <- data.frame(p.id = c(1, 1, 2, 2), period = c(1, 2, 1, 2), y = y, x.t = x.t,
                        z = c(z[1], z[1], z[2], z[2]))

  test.long <- WideToLong(wide.df, "p.id", "y", time.varying.bases = c("x.t"))
  checkTrue(all(test.long == long.df))
}

