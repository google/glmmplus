\name{CrossValidate}
\alias{CrossValidate}
\alias{crossValidate}
\alias{crossvalidate}
\alias{cross.validate}
\title{Cross Validation of gfo objects}
\description{ This is an elimination procedure for ...}

\usage{
CrossValidate(gfo, k.folds = 5, granularity = .005)
}

\arguments{
  \item{gfo}{
	  A gfo object from BackwardEliminate, ForwardSelect, or FitModel. 
  }
  \item{k.folds} { Number of cross validation splits }
  \item{granularity} { For binary case, how much to increment the threshold by}
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

}
\references{
	 Stef van Buuren, Karin Groothuis-Oudshoorn (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67. URL
	   http://www.jstatsoft.org/v45/i03/.
}
