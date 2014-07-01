\name{GetVariableImpacts}
\alias{GetVariableImpacts}
\alias{get.variable.impacts}
\alias{getVariableImpacts}
\title{Variable Impact Analysis after Cross Validation}
\description{Currently works with only a mids data object and family = binomial!}

\usage{
GetVariableImpacts(cv)
}

\arguments{
  \item{cv}{
	 A "gfo.cxv" object created from CrossValidate(). 
  }
}
\details{
The procedure works like this...
}
\examples{

data(testdata)
    
# creating a Muliply Imputed Data Set (mids) object
mids <- ImputeData(testdata, m = 5, maxit = 5)
gfo <- BackwardEliminate(y ~ x + w + z, data = mids)
cv <- CrossValidate(gfo)
impact <- GetVariableImpacts(cv)
impact
}
\references{
	 Stef van Buuren, Karin Groothuis-Oudshoorn (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67. URL
	   http://www.jstatsoft.org/v45/i03/.
}
