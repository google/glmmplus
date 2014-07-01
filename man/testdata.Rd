\name{testdata}
\alias{testdata}
\alias{TestData}
\title{A test data set for the glmmplus package}
\description{ Completely made up. y or y.binary are meant to be the response variables.
x, w, and z are continuous predictors, and factor.1 and factor.2 are categorical predictors.
The categorical predictors are also used as the id variable corresponding to a random effect
for testing of the lme4 based function calls.}

\usage{
data(testdata)
}

\details{
The data set was generated using R's random number generators.
}
\examples{
data(testdata)
head(testdata)
}

