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
advantage that they seamlessly work with the adjusted p-values from
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

After installing, load the glmmplus package. The package comes with the
`testdata` data frame that is useful for demonstrating package
functionality.

    library(glmmplus)
    # A sample data set with missing values
    head(testdata, 3)

    ##           x          z          w factor.1 factor.2         y y.binary
    ## 1 -1.207066  0.9847800 -1.2053334       15       30 -6.083643        0
    ## 2        NA -1.2247379  0.3014667        4       10  1.324600        1
    ## 3  1.084441  0.7097262 -1.5391452        2       15 -2.929349        0

Creating multiply imputed data sets
-----------------------------------

glmmplus provides a light wrapper around the mice package's `mice`
function that offers a single benefit: an argument called "droplist"
that allows the user to specify variables that will not be used in the
imputation models. This is because factor variables with many levels
often slow down the imputation process to unacceptable levels.

    my.mids <- ImputeData(testdata, m = 5, maxit = 5,
                          droplist = c("factor.1", "factor.2"))

The `complete` function has been exported from the mice package, and is
used below to create a data frame with no missing values based on a
single imputation.

    # a single imputation
    complete.df <- complete(my.mids, 1)
    head(complete.df, 3)

    ##            x          z          w factor.1 factor.2         y y.binary
    ## 1 -1.2070657  0.9847800 -1.2053334       15       30 -6.083643        0
    ## 2  0.4138689 -1.2247379  0.3014667        4       10  1.324600        1
    ## 3  1.0844412  0.7097262 -1.5391452        2       15 -2.929349        0

The model fitting functions in glmmplus wrap functions including
`base::glm`, `lme4::glmer`, and `nlme::lme` to provide a common
interface to generalized linear mixed modeling (including the cases
where there is no random factor). The output of these functions have
output class called "gfo," which stands for "generalized fitted object."

Forward selection with spline terms
-----------------------------------

Below is an example of a forward selection procedure on a data set with
no missing values. Note how the grouped term is associated with a single
p-value of the joint hypothesis.

    gfo.complete <- ForwardSelect(y.binary ~ ns(x, df = 2) + w + z,
                                  family = binomial, data = complete.df)
    summary(gfo.complete)

    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ w + ns(x, df = 2) 
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##             Term Wald p-value Unstandardized coefs Standardized coefs
    ## 1      Intercept         <NA>               -1.998                 NA
    ## 2              w       0.0000                0.811              0.761
    ## 3 ns(x, df = 2)1       0.0359                2.814                 NA
    ## 4 ns(x, df = 2)2   pval above                1.084                 NA
    ##   odds.ratio
    ## 1      0.136
    ## 2      2.250
    ## 3     16.676
    ## 4      2.956
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## ----Binomial response variable.----
    ## No well-defined R-squared is available.
    ## 
    ## McFadden's generalized R-squared: 0.1187 
    ## 
    ## Cox & Snell pseudo R-squared:  0.0964

Note that the family argument is currently set up to take functions and
not the usual string argument ([Issue
2](https://github.com/baogorek/glmmplus/issues/2)).

A mids object can be used in place of the data frame above, and all
inference is done on the Rubin's Rules adjusted p-values.

    gfo.mi <- ForwardSelect(y.binary ~ ns(x, df = 2) + w + z,
                                 family = binomial, data = my.mids)
    summary(gfo.mi)

    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ w 
    ## 
    ## Note: Multiple Imputation of missing values performed
    ## P-values adjusted according to Rubin's Rules
    ## 
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##        Term Wald p-value Unstandardized coefs Standardized coefs
    ## 1 Intercept         <NA>               -0.660                 NA
    ## 2         w       0.0000                0.818              0.767
    ##   odds.ratio
    ## 1      0.517
    ## 2      2.266
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## ----Binomial response variable.----
    ## No well-defined R-squared is available.
    ## 
    ## McFadden's generalized R-squared: 0.1063 
    ## 
    ## Cox & Snell pseudo R-squared:  0.0857

Incorporating random effects
----------------------------

Adding a random effect via lme4's syntax will automatically switch the
model fitting procedure to use `glmer`, and it works with mids objects
as well:

    # Backwards elimination for mixed (fixed and random) models
    ran.eff.gfo <- ForwardSelect(y.binary ~ (1 | factor.1) + x + w + z,
                                 family = binomial, data = my.mids)
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
    ##  factor.1 (Intercept) 1.51    
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##        Term Wald p-value Unstandardized coefs Standardized coefs
    ## 1 Intercept         <NA>               -0.913                 NA
    ## 2         w       0.0000                1.187              1.113
    ## 3         x       0.0444                0.374              0.389
    ##   odds.ratio
    ## 1      0.401
    ## 2      3.277
    ## 3      1.454
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## ----Binomial response variable.----
    ## No well-defined R-squared is available.
    ## 
    ## McFadden's generalized R-squared: 0.2921 
    ## 
    ## Cox & Snell pseudo R-squared:  0.2634

Incorporating within-subject time series structures
---------------------------------------------------

For a time-series model to be applied within the factor levels, add the
argument `ts.model = "ar1"` to the glmmplus model fitting functions. The
model fitting procedure then changes to `nlme::lme` and an
autoregressive structure will be applied within the single factor level
provided. In accordance with the `nlme::lme` function, there can only be
one random effects factor when using this option. Additionally, this is
only possible for gaussian response variables. time-series structures
are planned ([Issue 3](https://github.com/baogorek/glmmplus/issues/3)).

    ran.eff.ts.gfo <- ForwardSelect(y ~ (1 | factor.1) + x + w + z,
                                    family = gaussian, data = my.mids,
                                    ts.model = "ar1")
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
    ## 0.06 
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##        Term Wald p-value Unstandardized coefs Standardized coefs
    ## 1 Intercept         <NA>               -1.007                 NA
    ## 2         w       0.0000                1.532              0.466
    ## 3         x       0.0010                0.523              0.177
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## Random effects present. Using McFadden's R-Square with intercept-only null model
    ## R-squared: 0.5221

If you prefer backwards sequential selection, then use the
\`BackwardEliminate' function:

    ran.eff.ts.gfo.be <- BackwardEliminate(y ~ (1 | factor.1) + x + w + z,
                                           family = gaussian, data = my.mids,
                                           ts.model = "ar1")
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
    ## 0.06 
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##        Term Wald p-value Unstandardized coefs Standardized coefs
    ## 1 Intercept         <NA>               -1.007                 NA
    ## 2         x       0.0010                0.523              0.177
    ## 3         w       0.0000                1.532              0.466
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## Random effects present. Using McFadden's R-Square with intercept-only null model
    ## R-squared: 0.5221

Fitting a model without variable selection
------------------------------------------

Or, if you'd prefer not to use sequential variable selection methods at
all, use the `FitModel` function:

    ran.eff.ts.gfo.all <- FitModel(y ~ (1 | factor.1) + factor.2 + x + w + z,
                                   family = gaussian, data = my.mids,
                                   ts.model = "ar1")
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
    ## -0.04 
    ## 
    ## ### Type II analysis of fitted terms ###
    ## 
    ##          Term Wald p-value Unstandardized coefs Standardized coefs
    ## 1   Intercept         <NA>                0.410                 NA
    ## 2   factor.22       0.0000               -0.850                 NA
    ## 3   factor.23   pval above               -1.220                 NA
    ## 4   factor.24   pval above               -0.749                 NA
    ## 5   factor.25   pval above               -1.268                 NA
    ## 6   factor.26   pval above               -0.234                 NA
    ## 7   factor.27   pval above               -2.410                 NA
    ## 8   factor.28   pval above               -1.257                 NA
    ## 9   factor.29   pval above                2.595                 NA
    ## 10 factor.210   pval above               -2.624                 NA
    ## 11 factor.211   pval above               -1.011                 NA
    ## 12 factor.212   pval above               -3.821                 NA
    ## 13 factor.213   pval above               -1.095                 NA
    ## 14 factor.214   pval above               -0.239                 NA
    ## 15 factor.215   pval above               -4.654                 NA
    ## 16 factor.216   pval above               -2.693                 NA
    ## 17 factor.217   pval above                2.343                 NA
    ## 18 factor.218   pval above               -2.311                 NA
    ## 19 factor.219   pval above               -3.562                 NA
    ## 20 factor.220   pval above               -0.786                 NA
    ## 21 factor.221   pval above               -2.510                 NA
    ## 22 factor.222   pval above                0.805                 NA
    ## 23 factor.223   pval above                0.464                 NA
    ## 24 factor.224   pval above               -0.631                 NA
    ## 25 factor.225   pval above               -3.239                 NA
    ## 26 factor.226   pval above               -0.841                 NA
    ## 27 factor.227   pval above               -0.684                 NA
    ## 28 factor.228   pval above               -4.886                 NA
    ## 29 factor.229   pval above               -0.641                 NA
    ## 30 factor.230   pval above               -1.427                 NA
    ## 31          x       0.0000                0.606              0.204
    ## 32          w       0.0000                1.487              0.453
    ## 33          z       0.8496                0.012              0.004
    ## 
    ##  ### Fit Statistics ###
    ## 
    ## Random effects present. Using McFadden's R-Square with intercept-only null model
    ## R-squared: 0.8314

Notice that the example above also uses a fixed factor, where like the
spline term in example above, only a single p-value from the grouped
hypothesis test is presented.

References
==========

\[1\]: Rubin, D. B. (1987). Multiple Imputation for Nonresponse in
Surveys. New York: Wiley

\[2\]: Stef van Buuren, Karin Groothuis-Oudshoorn (2011). mice:
Multivariate Imputation by Chained Equations in R. Journal of
Statistical Software, 45(3), 1-67. URL
<http://www.jstatsoft.org/v45/i03/>.

\[3\]: Tibshirani, Robert (1996). Regression Shrinkage and Selection via
the Lasso. Journal of the Royal Statistical Society. Series B
(Methodological), Volume 58, Issue 1, 267-288.
