\name{MakeSplineStr}
\alias{MakeSplineStr}
\alias{MakeSplineString}
\alias{SplineString}
\alias{SplineStr}
\alias{MakeSpline}
\title{spline helper function}
\description{ This is  ...}

\usage{
MakeSplineStr(var.name, data, interior.knots)
}

\arguments{
  \item{var.name}{A variable name (character string)}
  \item{data}{a data frame}
  \item{interior.knots}{a numeric vector containing the interior knots of a 
	  parametric spline}
}
\details{
The procedure works like this...
}
\examples{
data(testdata)
spline.fn <- MakeSplineStr("x", testdata, interior.knots = c(-.5, .5))
}
