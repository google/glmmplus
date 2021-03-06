\name{ForwardSelect}
\alias{ForwardSelect}
\alias{forwardselect}
\alias{fs}
\title{Forward Selection for Generalized (Mixed) Linear Models with Missing Data}
\description{ This is an elimination procedure for ...}

\usage{
ForwardSelect(formula, data, cutoff, family)
}

\arguments{
  \item{formula}{
          A formula which may contain random effects according to the lme4 package's specification.      
  }                                                                                 
  \item{data}{
          Either a mids object from the mice package, or a data frame.
  }
  \item{cutoff}{
     The alpha level which determines the stopping rule. Once all remaining model terms fall below this value, 
     the procedure terminates.
  }
  \item{family}{
     Any family accepted by glm or lmer. Do not use quotation marks.
  }
}    

\details{
The procedure works like this...
}
\examples{
data(missing)

# A sample data set with missing values
head(missing)

# creating a Muliply Imputed Data Set (mids) object
mids <- ImputeData(missing, m = 5, maxit = 5)

# a single imputation
complete <- complete(mids)

# Backwards elimination for fixed effect models 
ForwardSelect(y ~ x + w + z, data = complete)
ForwardSelect(y ~ x + w + z, data = mids)

# Backwards elimination for mixed (fixed and random) models
ForwardSelect(y ~ (1 | factor.1) + x + w + z, data = complete)
ForwardSelect(y ~ (1 | factor.1) + x + w + z, data = mids)
}
\references{
Douglas Bates and Martin Maechler (2010). lme4: Linear mixed-effects models using S4 classes. R package
 version 0.999375-37. http://CRAN.R-project.org/package=lme4

Stef van Buuren, Karin Groothuis-Oudshoorn (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67. URL http://www.jstatsoft.org/v45/i03/.
}
