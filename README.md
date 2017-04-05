glmmplus: generalized linear mixed models with multiple imputation
==================================================================

[![Build
Status](https://travis-ci.org/baogorek/glmmplus.svg?branch=master)](https://travis-ci.org/glmmplus)

[![Coverage
Status](https://img.shields.io/codecov/c/github/baogorek/glmmplus/master.svg)](https://codecov.io/github/baogorek/glmmplus?branch=master)

Overview
--------

glmmplus was based on the notion that missing data is a fundamental part
of nearly any analysis and that the flagship procedures that allow
fitting of linear models and variable selection should include better
options than ejecting rows for dealing with missing values. Due to the
availability and general acceptance of multiple imputation via the mice
package, glmnet provides an interface that takes mids (multiply imputed
data sets) and data frame objects interchangeably as part of its
interface.

Though having fallen out of favor since the arrival of LASSO, sequential
variable selection methods that are based on p-values have the advantage
that they seemlessly work with the adjusted p-values from multiple
imputation (Rubin's Rules). Thus, a large part of glmmplus's code is
creating sequential variable selection routines based on p-values.

This package was developed while I was employed at Google, so it follows
the Google style guide. An older version is still up at the github repo
google/glmmplus.

Installation
------------

An easy way to install glmmplus is via the devtools package:

    devtools::install_github("baogorek/glmmplus")

Creating multiply imputed data sets
-----------------------------------

glmmplus provides a light wrapper around the mice function (from the
mice) package that offers one main benefit: an argument called
"droplist" that allows the user to specify variables that will not be
used in the imputation models. (Factor variables with many levels often
slow down the imputation process to unacceptable levels.)

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
    ## 2 -1.1719483 -1.2247379  0.3014667        4       10  1.3246005        1
    ## 3  1.0844412  0.7097262 -1.5391452        2       15 -2.9293494        0
    ## 4 -2.3456977 -1.4191629  0.6353707        4       25  1.1263308        1
    ## 5  0.4291247  1.7826079  0.7029518        4        7  1.8670865        1
    ## 6  0.5060559 -0.2434447 -1.9058829        3        6 -0.9476595        0

The model fitting functions in glmmplus are wrap functions like
base::glm, lme4::glmer, and nlme::lme to provide a common interface. The
output of these functions have output of class "gfo" ("generalized
fitted object").

    # Backwards elimination for fixed effect models
    gfo.complete <- ForwardSelect(y.binary ~ ns(x, df = 2) + w + z,
                                  family = binomial, data = complete.df)

    ## 
    ##   **** Stepwise Iteration 1 ****
    ## Testing individual terms...
    ## --------
    ## var:  ns(x, df = 2) , pvalue:  0.0002  
    ## var:  w , pvalue:  0.0000  
    ## var:  z , pvalue:  0.3684  
    ## Forward selection iteration 1 ,adding term w to model
    ## 
    ## 
    ## Iteration 1 formula:
    ##  y.binary ~ w 
    ## 
    ##   **** Stepwise Iteration 2 ****
    ## Testing individual terms...
    ## --------
    ## var:  ns(x, df = 2) , pvalue:  0.0008  
    ## var:  z , pvalue:  0.2060  
    ## Forward selection iteration 2 ,adding term ns(x, df = 2) to model
    ## 
    ## 
    ## Iteration 2 formula:
    ##  y.binary ~ w + ns(x, df = 2) 
    ## 
    ##   **** Stepwise Iteration 3 ****
    ## Testing individual terms...
    ## --------
    ## var:  z , pvalue:  0.5068  
    ## No further terms have p-values < 0.05 
    ## 
    ## 
    ## Iteration 3 formula:
    ##  y.binary ~ w + ns(x, df = 2) 
    ## Generalized Fitted Object 'gfo'
    ## 
    ## R Formula Syntax: y.binary ~ w + ns(x, df = 2)

    class(gfo.complete)

    ## [1] "gfo"

    coef(gfo.complete)

    ##    (Intercept)              w ns(x, df = 2)1 ns(x, df = 2)2 
    ##     -2.3890725      0.8103394      3.7054872      1.7228632

Note that the family argument is currently set up to take functions and
not the usual string argment (TODO: change that). Also note that
throughout the forward selection process, the spline terms were
considered together. This is possible with LASSO but is typically not
offered in popular implementations (e.g., glmnet).

As discussed in the overview, a mids object can be used in place of a
data frame, and all inference is done on the adjusted p-values from
application of "Rubin's Rules."

    gfo.missing <- ForwardSelect(y.binary ~ ns(x, df = 2) + w + z,
                                 family = binomial, data = my.mids)

    ## 
    ##   **** Stepwise Iteration 1 ****
    ## Testing individual terms...
    ## --------
    ## var:  ns(x, df = 2) , pvalue:  0.0283  
    ## var:  w , pvalue:  0.0000  
    ## var:  z , pvalue:  0.7692  
    ## Forward selection iteration 1 ,adding term w to model
    ## 
    ## 
    ## Iteration 1 formula:
    ##  y.binary ~ w 
    ## 
    ##   **** Stepwise Iteration 2 ****
    ## Testing individual terms...
    ## --------
    ## var:  ns(x, df = 2) , pvalue:  0.0446  
    ## var:  z , pvalue:  0.5705  
    ## Forward selection iteration 2 ,adding term ns(x, df = 2) to model
    ## 
    ## 
    ## Iteration 2 formula:
    ##  y.binary ~ w + ns(x, df = 2) 
    ## 
    ##   **** Stepwise Iteration 3 ****
    ## Testing individual terms...
    ## --------
    ## var:  z , pvalue:  0.7594  
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

    class(gfo.missing)

    ## [1] "gfo"

    coef(gfo.missing)

    ##    (Intercept)              w ns(x, df = 2)1 ns(x, df = 2)2 
    ##     -1.5450807      0.8135906      2.0873039      1.6640041

Note that after the imputation, the spline term collectively fell below
the default threshold of .05 and was not included in the final model.

Incorporating random effects
----------------------------

Adding a random effect via lme4's syntax will automatically switch the
model fitting procedure to glmer, and it works with mids objects as
well:

    # Backwards elimination for mixed (fixed and random) models
    ran.eff.gfo <- ForwardSelect(y.binary ~ (1 | factor.1) + x + w + z,
                                 family = binomial, data = my.mids)

    ## 
    ##   **** Stepwise Iteration 1 ****
    ## Testing individual terms...
    ## --------
    ## var:  x , pvalue:  0.0051  
    ## var:  w , pvalue:  0.0000  
    ## var:  z , pvalue:  0.9113  
    ## Forward selection iteration 1 ,adding term w to model
    ## 
    ## 
    ## Iteration 1 formula:
    ##  y.binary ~ (1 | factor.1) + w 
    ## 
    ##   **** Stepwise Iteration 2 ****
    ## Testing individual terms...
    ## --------
    ## var:  x , pvalue:  0.0143  
    ## var:  z , pvalue:  0.5435  
    ## Forward selection iteration 2 ,adding term x to model
    ## 
    ## 
    ## Iteration 2 formula:
    ##  y.binary ~ (1 | factor.1) + w + x 
    ## 
    ##   **** Stepwise Iteration 3 ****
    ## Testing individual terms...
    ## --------
    ## var:  z , pvalue:  0.6564  
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

    coef(ran.eff.gfo)

    ## (Intercept)           w           x 
    ##  -0.9272539   1.2035952   0.4541231 
    ## attr(,"class")
    ## [1] "coef.mer"
