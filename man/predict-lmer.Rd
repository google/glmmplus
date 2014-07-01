\name{predict.mer}
\alias{predict.mer}
\alias{predict.lmer}
\title{Prediction function for Generalized (Mixed) Linear Models with Missing Data}
\description{Predicts given an lmer object}

\usage{
predict(obj, newdata)
}

\arguments{
  \item{obj}{
       An lmer object.
  }
  \item{newdata} {either a data frame or a mids object from the mice package}
}
\details{
The procedure works like this...
}
\examples{

  model <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
  preds <- predict(model, se.fit = TRUE, newdata = sleepstudy)
  plot(preds, plot.smean = TRUE, main = "Gaussian GLMM predictions")

  # http://www.statistik.lmu.de/~greven/alda/examples_mixed_models.R
  I <- 100
  n.i <- 8
  subject <- rep(1:I,each = n.i)
  random.intercept <- rep(rnorm(I),each = n.i)
  time <- rep(1:n.i, I)
  beta.0 <- beta.1 <- 1
  lambda <- exp(beta.0 + random.intercept + beta.1 * time)
  Y <- rpois(I * n.i,lambda)
  data <- data.frame(Y, time, subject)

  model <- glmer(Y ~ time + (1 | subject), data, family = poisson)
  preds <- predict(model, se.fit = TRUE, newdata = data)
  plot(preds, main = "Poisson GLMM predictions")

  # http://www.ats.ucla.edu/stat/r/dae/melogit.htm
  hdp <- read.csv("http://www.ats.ucla.edu/stat/data/hdp.csv")
  model <-
    glmer(remission ~ IL6 + CRP + CancerStage + LengthofStay + Experience +
          (1 | DID), data = hdp, family = binomial, nAGQ = 10)
  preds <- predict(model, se.fit = TRUE, newdata = hdp)
  plot(preds ,  plot.smean = TRUE,
       main = "Binomial GLMM predictions")
}

\references{
         Stef van Buuren, Karin Groothuis-Oudshoorn (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67. URL
           http://www.jstatsoft.org/v45/i03/.

  TODO: finish citation. http://www.statistik.lmu.de/~greven/alda/examples_mixed_models.R

  TODO: finish citaiton. http://www.ats.ucla.edu/stat/r/dae/melogit.htm

}
