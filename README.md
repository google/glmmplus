glmmplus: generalized linear mixed models, model selection, and multiple imputation
===================================================================================

[![Build
Status](https://travis-ci.org/baogorek/glmmplus.svg?branch=master)](https://travis-ci.org/glmmplus)

[![Coverage
Status](https://img.shields.io/codecov/c/github/baogorek/glmmplus/master.svg)](https://codecov.io/github/baogorek/glmmplus?branch=master)

Overview
--------

Missing data is a key feature of nearly any modern regression analysis,
but many of R's flagship modeling tools default to ejecting rows. The
glmmplus package uses Wald p-values for estimation and model selection
in conjunction with the lme4, nlme, and stats packages to enable
consistent application of Multiple imputation (\[Rubin, 1983\]\[1\]), a
generally accepted and principled method of incorporating additional
uncertainty into an analysis with imputation. As the mice package (\[van
Buuren and Groothuis-Oudshoorn, 2011\]\[2\]) is a popular and trusted
implementation of multiple imputation in R, glmmplus builds upon the
mids (multiply imputed data set) interface, allowing mids objects and
data frames to be used interchangeably.

The glmmplus package provides sequential model selection of fixed
effects using Rubin's adjusted p-values. Though having fallen out of
favor since the popularization of LASSO (\[Tibshirani 1996\]\[3\]),
sequential variable selection methods based on p-values have the
advantage that they seemlessly work with the adjusted p-values from
multiple imputation, and easily allow for grouped terms and hierarchies
to be considered. Thus, a large part of glmmplus's code is creating
sequential variable selection routines based on p-values.

This package was originally developed within Google's People Analytics
department, so the output has been designed to satisfy a social
scientist's expectations. Additionally, it follows the Google style
guide. While the original version is still up at Google's github repo
(google/glmmplus), there is no one on the other side to take pull
requests and this version is dated, full of bugs, and should not be
used.

Installation
------------

An easy way to install glmmplus is via the devtools package:

    devtools::install_github("baogorek/glmmplus")

Creating multiply imputed data sets
-----------------------------------

glmmplus provides a light wrapper around the mice package's `mice`
function that offers a single benefit: an argument called "droplist"
that allows the user to specify variables that will not be used in the
imputation models. This is becuase factor variables with many levels
often slow down the imputation process to unacceptable levels.

    library(glmmplus)
    # A sample data set with missing values
    head(testdata)

    ##            x          z          w factor.1 factor.2          y y.binary
    ## 1 -1.2070657  0.9847800 -1.2053334       15       30 -6.0836429        0
    ## 2         NA -1.2247379  0.3014667        4       10  1.3246005        1
    ## 3  1.0844412  0.7097262 -1.5391452        2       15 -2.9293494        0
    ## 4 -2.3456977         NA  0.6353707        4       25  1.1263308        1
    ## 5  0.4291247  1.7826079  0.7029518        4        7  1.8670865        1
    ## 6  0.5060559 -0.2434447 -1.9058829        3        6 -0.9476595        0

    # creating a Muliply Imputed Data Set (mids) object
    my.mids <- ImputeData(testdata, m = 5, maxit = 5,
                          droplist = c("factor.1", "factor.2"))

    ## 
    ##  iter imp variable
    ##   1   1  x  z
    ##   1   2  x  z
    ##   1   3  x  z
    ##   1   4  x  z
    ##   1   5  x  z
    ##   2   1  x  z
    ##   2   2  x  z
    ##   2   3  x  z
    ##   2   4  x  z
    ##   2   5  x  z
    ##   3   1  x  z
    ##   3   2  x  z
    ##   3   3  x  z
    ##   3   4  x  z
    ##   3   5  x  z
    ##   4   1  x  z
    ##   4   2  x  z
    ##   4   3  x  z
    ##   4   4  x  z
    ##   4   5  x  z
    ##   5   1  x  z
    ##   5   2  x  z
    ##   5   3  x  z
    ##   5   4  x  z
    ##   5   5  x  z

    # a single imputation
    complete.df <- complete(my.mids)
    head(complete.df)

    ##            x          z          w factor.1 factor.2          y y.binary
    ## 1 -1.2070657  0.9847800 -1.2053334       15       30 -6.0836429        0
    ## 2 -0.9532787 -1.2247379  0.3014667        4       10  1.3246005        1
    ## 3  1.0844412  0.7097262 -1.5391452        2       15 -2.9293494        0
    ## 4 -2.3456977 -1.7579991  0.6353707        4       25  1.1263308        1
    ## 5  0.4291247  1.7826079  0.7029518        4        7  1.8670865        1
    ## 6  0.5060559 -0.2434447 -1.9058829        3        6 -0.9476595        0

The model fitting functions in glmmplus wrap functions including
`base::glm`, `lme4::glmer`, and `nlme::lme` to provide a common
interface to the provided functionality. The output of these functions
have output class called "gfo," which stands for "generalized fitted
object."

Below is an example of a forward selection procedure on a data set with
no missing values. Note how the grouped term is associated with a single
p-value of the joint hypothesis.

    gfo.complete <- ForwardSelect(y.binary ~ ns(x, df = 2) + w + z,
                                  family = binomial, data = complete.df)

    ## 
    ##   **** Stepwise Iteration 1 ****
    ## Testing individual terms...
    ## --------
    ## var:  ns(x, df = 2) , pvalue:  0.0146  
    ## var:  w , pvalue:  0.0000  
    ## var:  z , pvalue:  0.3134  
    ## Forward selection iteration 1 ,adding term w to model
    ## 
    ## 
    ## Iteration 1 formula:
    ##  y.binary ~ w 
    ## 
    ##   **** Stepwise Iteration 2 ****
    ## Testing individual terms...
    ## --------
    ## var:  ns(x, df = 2) , pvalue:  0.0230  
    ## var:  z , pvalue:  0.1048  
    ## Forward selection iteration 2 ,adding term ns(x, df = 2) to model
    ## 
    ## 
    ## Iteration 2 formula:
    ##  y.binary ~ w + ns(x, df = 2) 
    ## 
    ##   **** Stepwise Iteration 3 ****
    ## Testing individual terms...
    ## --------
    ## var:  z , pvalue:  0.2095  
    ## No further terms have p-values < 0.05 
    ## 
    ## 
    ## Iteration 3 formula:
    ##  y.binary ~ w + ns(x, df = 2) 
    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ w + ns(x, df = 2) 
    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ w + ns(x, df = 2)

    summary(gfo.complete)

    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ w + ns(x, df = 2) 
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##             Term Wald p-value Unstandardized coefs Standarized coefs
    ## 1      Intercept         <NA>               -1.279                NA
    ## 2              w       0.0000                0.820             0.769
    ## 3 ns(x, df = 2)1       0.0230                1.500                NA
    ## 4 ns(x, df = 2)2   pval above                1.356                NA
    ##   odds.ratio
    ## 1      0.278
    ## 2      2.270
    ## 3      4.482
    ## 4      3.881
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## ----Binomial response varible.----
    ## No well-defined R-squared is available.
    ## 
    ## McFadden's generalized R-squared: 0.1198 
    ## 
    ## Cox & Snell pseudo R-squared:  0.0973

Note that the family argument is currently set up to take functions and
not the usual string argment ([Issue
2](https://github.com/baogorek/glmmplus/issues/2)).

A mids object can be used in place of the data frame above, and all
inference is done on the Rubin's Rules adjusted p-values.

    gfo.mi <- ForwardSelect(y.binary ~ ns(x, df = 2) + w + z,
                                 family = binomial, data = my.mids)

    ## 
    ##   **** Stepwise Iteration 1 ****
    ## Testing individual terms...
    ## --------
    ## var:  ns(x, df = 2) , pvalue:  0.0097  
    ## var:  w , pvalue:  0.0000  
    ## var:  z , pvalue:  0.9765  
    ## Forward selection iteration 1 ,adding term w to model
    ## 
    ## 
    ## Iteration 1 formula:
    ##  y.binary ~ w 
    ## 
    ##   **** Stepwise Iteration 2 ****
    ## Testing individual terms...
    ## --------
    ## var:  ns(x, df = 2) , pvalue:  0.0158  
    ## var:  z , pvalue:  0.8523  
    ## Forward selection iteration 2 ,adding term ns(x, df = 2) to model
    ## 
    ## 
    ## Iteration 2 formula:
    ##  y.binary ~ w + ns(x, df = 2) 
    ## 
    ##   **** Stepwise Iteration 3 ****
    ## Testing individual terms...
    ## --------
    ## var:  z , pvalue:  0.9738  
    ## No further terms have p-values < 0.05 
    ## 
    ## 
    ## Iteration 3 formula:
    ##  y.binary ~ w + ns(x, df = 2) 
    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ w + ns(x, df = 2) 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules
    ## 
    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ w + ns(x, df = 2) 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules

    summary(gfo.mi)

    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ w + ns(x, df = 2) 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules
    ## 
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##             Term Wald p-value Unstandardized coefs Standarized coefs
    ## 1      Intercept         <NA>               -1.354                NA
    ## 2              w       0.0000                0.816             0.765
    ## 3 ns(x, df = 2)1       0.0158                1.691                NA
    ## 4 ns(x, df = 2)2   pval above                1.533                NA
    ##   odds.ratio
    ## 1      0.258
    ## 2      2.261
    ## 3      5.425
    ## 4      4.632
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## ----Binomial response varible.----
    ## No well-defined R-squared is available.
    ## 
    ## McFadden's generalized R-squared: 0.1222 
    ## 
    ## Cox & Snell pseudo R-squared:  0.0994

Incorporating random effects
----------------------------

Adding a random effect via lme4's syntax will automatically switch the
model fitting procedure to use `glmer`, and it works with mids objects
as well:

    # Backwards elimination for mixed (fixed and random) models
    ran.eff.gfo <- ForwardSelect(y.binary ~ (1 | factor.1) + x + w + z,
                                 family = binomial, data = my.mids)

    ## 
    ##   **** Stepwise Iteration 1 ****
    ## Testing individual terms...
    ## --------
    ## var:  x , pvalue:  0.0009  
    ## var:  w , pvalue:  0.0000  
    ## var:  z , pvalue:  0.8563  
    ## Forward selection iteration 1 ,adding term w to model
    ## 
    ## 
    ## Iteration 1 formula:
    ##  y.binary ~ (1 | factor.1) + w 
    ## 
    ##   **** Stepwise Iteration 2 ****
    ## Testing individual terms...
    ## --------
    ## var:  x , pvalue:  0.0013  
    ## var:  z , pvalue:  0.8580  
    ## Forward selection iteration 2 ,adding term x to model
    ## 
    ## 
    ## Iteration 2 formula:
    ##  y.binary ~ (1 | factor.1) + w + x 
    ## 
    ##   **** Stepwise Iteration 3 ****
    ## Testing individual terms...
    ## --------
    ## var:  z , pvalue:  0.9701  
    ## No further terms have p-values < 0.05 
    ## 
    ## 
    ## Iteration 3 formula:
    ##  y.binary ~ (1 | factor.1) + w + x 
    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ (1 | factor.1) + w + x 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules
    ## 
    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ (1 | factor.1) + w + x 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules

    summary(ran.eff.gfo)

    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ (1 | factor.1) + w + x 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules
    ## 
    ## 
    ## 
    ## ###  Random Effects  ###
    ## 
    ##  Groups   Name        Std.Dev.
    ##  factor.1 (Intercept) 1.52    
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##        Term Wald p-value Unstandardized coefs Standarized coefs odds.ratio
    ## 1 Intercept         <NA>               -0.921                NA      0.398
    ## 2         w       0.0000                1.203             1.128      3.330
    ## 3         x       0.0013                0.419             0.435      1.520
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## ----Binomial response varible.----
    ## No well-defined R-squared is available.
    ## 
    ## McFadden's generalized R-squared: 0.2961 
    ## 
    ## Cox & Snell pseudo R-squared:  0.2678

For a time-series model to be applied within the factor levels, add the
argument `ts.model = "ar1"` to the glmmplus model fitting functions. The
model fitting procedure then changes to `nlme::lme` and an
autoregressive structure will be applied within the single factor level
provided. In accordance with the `nlme::lme1` function, there can only
be one random effects factor when using this option. Additionally, this
is only possible for gaussian response variables. time-series structures
are planned ([Issue 3](https://github.com/baogorek/glmmplus/issues/3)).

    ran.eff.ts.gfo <- ForwardSelect(y ~ (1 | factor.1) + x + w + z,
                                    family = gaussian, data = my.mids,
                                    ts.model = "ar1")

    ## 
    ##   **** Stepwise Iteration 1 ****
    ## Testing individual terms...
    ## --------
    ## var:  x , pvalue:  0.0000  
    ## var:  w , pvalue:  0.0000  
    ## var:  z , pvalue:  0.6709  
    ## Forward selection iteration 1 ,adding term w to model
    ## 
    ## 
    ## Iteration 1 formula:
    ##  y ~ (1 | factor.1) + w 
    ## 
    ##   **** Stepwise Iteration 2 ****
    ## Testing individual terms...
    ## --------
    ## var:  x , pvalue:  0.0000  
    ## var:  z , pvalue:  0.4288  
    ## Forward selection iteration 2 ,adding term x to model
    ## 
    ## 
    ## Iteration 2 formula:
    ##  y ~ (1 | factor.1) + w + x 
    ## 
    ##   **** Stepwise Iteration 3 ****
    ## Testing individual terms...
    ## --------
    ## var:  z , pvalue:  0.6600  
    ## No further terms have p-values < 0.05 
    ## 
    ## 
    ## Iteration 3 formula:
    ##  y ~ (1 | factor.1) + w + x 
    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y ~ (1 | factor.1) + w + x 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules
    ## 
    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y ~ (1 | factor.1) + w + x 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules

    summary(ran.eff.ts.gfo)

    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y ~ (1 | factor.1) + w + x 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules
    ## 
    ## 
    ## 
    ## ###  Random Effects  ###
    ## 
    ##     Groups        Name Variance
    ## 1  Subject (Intercept)      1.8
    ## 2 Residual                  2.0
    ## 
    ## AR1 parameter estimate from lme:
    ## 0.05 
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##        Term Wald p-value Unstandardized coefs Standarized coefs
    ## 1 Intercept         <NA>               -1.011                NA
    ## 2         w       0.0000                1.545             0.470
    ## 3         x       0.0000                0.515             0.174
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## Random effects present. Using McFadden's R-Square with intercept-only null model
    ## R-squared: 0.5226

If you prefer backwards sequential selection, then use the
\`BackwardEliminate' function:

    ran.eff.ts.gfo.be <- BackwardEliminate(y ~ (1 | factor.1) + x + w + z,
                                           family = gaussian, data = my.mids,
                                           ts.model = "ar1")

    ## 
    ##   **** Stepwise Iteration 1 ****
    ## Testing individual terms...
    ## --------
    ## var:  x , pvalue:  0.0000  
    ## var:  w , pvalue:  0.0000  
    ## var:  z , pvalue:  0.6600  
    ## Backward Elimination iteration 1 , dropping term z from model
    ## 
    ## 
    ## Iteration 1 formula:
    ##  y ~ (1 | factor.1) + x + w 
    ## 
    ##   **** Stepwise Iteration 2 ****
    ## Testing individual terms...
    ## --------
    ## var:  x , pvalue:  0.0000  
    ## var:  w , pvalue:  0.0000  
    ## All remaining terms have p-values < 0.05 
    ## 
    ## 
    ## Iteration 2 formula:
    ##  y ~ (1 | factor.1) + x + w 
    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y ~ (1 | factor.1) + x + w 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules

    summary(ran.eff.ts.gfo.be)

    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y ~ (1 | factor.1) + x + w 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules
    ## 
    ## 
    ## 
    ## ###  Random Effects  ###
    ## 
    ##     Groups        Name Variance
    ## 1  Subject (Intercept)      1.8
    ## 2 Residual                  2.0
    ## 
    ## AR1 parameter estimate from lme:
    ## 0.05 
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##        Term Wald p-value Unstandardized coefs Standarized coefs
    ## 1 Intercept         <NA>               -1.011                NA
    ## 2         x       0.0000                0.515             0.174
    ## 3         w       0.0000                1.545             0.470
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## Random effects present. Using McFadden's R-Square with intercept-only null model
    ## R-squared: 0.5226

Or, if you'd prefer not to use sequential variable selection methods at
all, use the `FitModel` function:

    ran.eff.ts.gfo.all <- FitModel(y ~ (1 | factor.1) + factor.2 + x + w + z,
                                   family = gaussian, data = my.mids,
                                   ts.model = "ar1")

    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y ~ (1 | factor.1) + factor.2 + x + w + z 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules

    summary(ran.eff.ts.gfo.all)

    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y ~ (1 | factor.1) + factor.2 + x + w + z 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules
    ## 
    ## 
    ## 
    ## ###  Random Effects  ###
    ## 
    ##     Groups        Name Variance
    ## 1  Subject (Intercept)      1.9
    ## 2 Residual                  1.2
    ## 
    ## AR1 parameter estimate from lme:
    ## -0.05 
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##          Term Wald p-value Unstandardized coefs Standarized coefs
    ## 1   Intercept         <NA>                0.388                NA
    ## 2   factor.22       0.0000               -0.790                NA
    ## 3   factor.23   pval above               -1.188                NA
    ## 4   factor.24   pval above               -0.753                NA
    ## 5   factor.25   pval above               -1.332                NA
    ## 6   factor.26   pval above               -0.239                NA
    ## 7   factor.27   pval above               -2.464                NA
    ## 8   factor.28   pval above               -1.345                NA
    ## 9   factor.29   pval above                2.650                NA
    ## 10 factor.210   pval above               -2.572                NA
    ## 11 factor.211   pval above               -0.959                NA
    ## 12 factor.212   pval above               -3.783                NA
    ## 13 factor.213   pval above               -1.101                NA
    ## 14 factor.214   pval above               -0.198                NA
    ## 15 factor.215   pval above               -4.556                NA
    ## 16 factor.216   pval above               -2.725                NA
    ## 17 factor.217   pval above                2.470                NA
    ## 18 factor.218   pval above               -2.288                NA
    ## 19 factor.219   pval above               -3.515                NA
    ## 20 factor.220   pval above               -0.821                NA
    ## 21 factor.221   pval above               -2.488                NA
    ## 22 factor.222   pval above                0.745                NA
    ## 23 factor.223   pval above                0.502                NA
    ## 24 factor.224   pval above               -0.601                NA
    ## 25 factor.225   pval above               -3.223                NA
    ## 26 factor.226   pval above               -0.847                NA
    ## 27 factor.227   pval above               -0.727                NA
    ## 28 factor.228   pval above               -4.852                NA
    ## 29 factor.229   pval above               -0.650                NA
    ## 30 factor.230   pval above               -1.380                NA
    ## 31          x       0.0000                0.581             0.196
    ## 32          w       0.0000                1.503             0.457
    ## 33          z       0.9908                0.001             0.000
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## Random effects present. Using McFadden's R-Square with intercept-only null model
    ## R-squared: 0.8304

Notice that the example above also uses a fixed factor, where like the
spline term in example above, only a single p-value from the grouped
hypothesis test is presented.

\[1\]: Rubin, D. B. (1987). Multiple Imputation for Nonresponse in
Surveys. New York: Wiley

\[2\]: Stef van Buuren, Karin Groothuis-Oudshoorn (2011). mice:
Multivariate Imputation by Chained Equations in R. Journal of
Statistical Software, 45(3), 1-67. URL
<http://www.jstatsoft.org/v45/i03/>.

\[3\]: Tibshirani, Robert (1996). Regression Shrinkage and Selection via
the Lasso. Journal of the Royal Statistical Society. Series B
(Methodological), Volume 58, Issue 1, 267-288.
