\name{GetVariableGroupImpact}
\alias{GetVariableGroupImpact}
\alias{GetVariableGroupImpacts}

\title{Model Implied Group Impact Analysis}
\description{Shows how the predicted values of a model change as a set of variables 
             simultaneously change in value.}

\usage{
GetVariableGroupImpact(gfo.obj, var.vec, low.val.list, high.val.list)
}

\arguments{
  \item{gfo.obj}{
	  A gfo object created by FitModel, BackwardEliminate, or Forward Select
  }
  \item{var.vec}{
	  A character vector of variable names
  }
  \item{low.val.list}{
          A list mapping each variable name (from var.vec) to its "low" value
  }
  \item{high.val.list}{
	  A list mapping each variable name (from var.vec) to its "high" value
  }
}
\details{
The procedure works like this...
}

\examples{
data(testdata)

# A sample data set with testdata values
head(testdata)

# Backwards elimination for fixed effect models 
gfo <- FitModel(y ~ x + w + z, data = testdata)
GetVariableGroupImpact(gfo, var.vec = c("w", "x"),
		       low.val.list = list(x = -1, w = -1),
		       high.val.list = list(x = 1, w = 1))


}
