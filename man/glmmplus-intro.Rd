\name{glmmplus}
\alias{glmmplus}
\title{glmmplus - An introduction}
\description{
Welcome to the glmmplus package. For documentation, please see
\url{go/glmmplus}.

Briefly, glmmplus is a composite R package which ties together functionality from lme4, mice, and the R stats functionality. 
A multiple imputation package, mice offers the computational strategy laid out in ... As of version 1.0, tests using mice did not 
work with mer class models created with the lme4 package, and only single tests were supported. We have extended the Wald testing
strategy of the mice package to work with a wider variety of models and variable selection algorithms.

Selection of random effects...

This package comes with a set of unit tests (see the \file{tests}).

The main functions of this package are FitModel, BackwardEliminate, and ForwardSelect.
For examples, please see the help pages corresponding to these functions.
}
