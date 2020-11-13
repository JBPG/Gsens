Gsens
================
Jean-Baptiste Pingault & Frank Dudbridge
20 August 2020

### Overview

#### [1: Introduction](#link1)

#### [2: One polygenic score case](#link2)

#### [3: Two polygenic scores](#link3)

# Introduction<a name="link1"></a>

The following script provides an example on how to run Gsens, a
genetically informed sensitivity analysis, in R. Please cite the
following papers where the concept was first proposed and then
implemented:

1.  Pingault, J.-B., O’Reilly, P. F., Schoeler, T., Ploubidis, G. B.,
    Rijsdijk, F., & Dudbridge, F. (2018). Using genetic data to
    strengthen causal inference in observational research. Nature
    Reviews Genetics, 19(9), 566–580.
    <https://doi.org/10.1038/s41576-018-0020-3>

2.  (UPDATE) Pingault, J.-B., Rijsdijk, F., Schoeler, T., Choi, S. W.,
    Kraphol, E., O’Reilly, P. F., & Dudbridge, F. (In Prep). Estimating
    the sensitivity of associations between risk factors and outcomes to
    shared genetic effects. BioRxiv.

In this tutorial, we illustrate the use of `Gsens` to estimate the role
of genetic confounding in explaining the associations between maternal
educational and three developmental outcomes in the child

  - educational achievement
  - BMI
  - ADHD

The following examples are based on the correlation matrix between
polygenic scores and variables, as shown below (also available in the
supplementary material of article 2).

**Correlations between study
variables**

|                        | Maternal Education | GCSE     | BMI      | ADHD     | EDU PS   | BMI PS   | ADHD PS  |
| ---------------------- | ------------------ | -------- | -------- | -------- | -------- | -------- | -------- |
| **Maternal Education** | 1                  | 0.3975   | \-0.0894 | \-0.1240 | 0.2909   | \-0.0839 | \-0.0630 |
| **GCSE**               | 0.3975             | 1        | \-0.0894 | \-0.3398 | 0.3462   | \-0.0885 | \-0.1103 |
| **BMI**                | \-0.0894           | \-0.0894 | 1        | \-0.0088 | \-0.0269 | 0.2524   | 0.0441   |
| **ADHD**               | \-0.1240           | \-0.3398 | \-0.0088 | 1        | \-0.0902 | 0.0820   | 0.1188   |
| **EDU PS**             | 0.2909             | 0.3462   | \-0.0269 | \-0.0902 | 1        | \-0.1856 | \-0.1865 |
| **BMI PS**             | \-0.0839           | \-0.0885 | 0.2524   | 0.0820   | \-0.1856 | 1        | 0.1413   |
| **ADHD PS**            | \-0.0630           | \-0.1103 | 0.0441   | 0.1188   | \-0.1865 | 0.1413   | 1        |

</br>

# Installation of Gsens

</br>

To begin, please install and load the Gsens package with the code below:

``` r
# Install devtools
install.packages("devtools")
library(devtools)

# Install Gsens
install_github("JBPG/Gsens")
library(Gsens)
```

Three functions are available:

`gsensY()`: sensivity analysis based on one polygenic score for the
outcome (Y)

`gsensX()`: sensivity analysis based on one polygenic score for the
exposure (X)

`gsensXY()`: sensivity analysis based on two polygenic scores for X and
Y

# One polygenic score case<a name="link2"></a>

## Observed scenario

In the following example, we will test the association between based on
maternal years of education (X) and child GCSE scores (Y) after
controlling for the best fitting polygenic score for years of education
estimated in the child.

To do so, we will use the “gsensY” function.

For this, we need to specify 5 parameters:

  - **rxy** = the observed phenotypic correlations between study
    variables (here: correlation between maternal years of education and
    child GCSE);

  - **rgx** = the observed correlations between phenotype X (maternal
    years of education) and the observed polygenic score for education
    (note that these should be partial correlations adjusting for sex,
    age, and principal components);

  - **rgy** = the observed correlations between phenotype Y (child GCSE)
    and the observed polygenic score for education (adjusted for sex,
    age, and principal components);

  - **n** = sample size;

  - **h2** = is the variance explained in the outcome, here by the
    observed polygenic score (hence why h2 is **rgy^2**).

Note that rxz, rgx, and rgy should be partial correlations adjusting for
sex, age, and principal components

``` r
gsensY(rxy=0.3975,
             rgx = 0.2909,
             rgy = 0.3462,
             n=3785,
             h2=0.3462^2)
```

    ##                       est    se      z      pvalue ci.lower ci.upper
    ## Adjusted Bxy        0.324 0.017 19.396  8.4077e-84    0.291    0.357
    ## Genetic confounding 0.073 0.003 22.881 7.1055e-116    0.067    0.080
    ## Total effect        0.397 0.015 26.479 1.6783e-154    0.368    0.427

The output provides:

  - **Adjusted Bxy** This is the standardized estimate of the
    relationship between X and Y, adjusted for G (i.e. the residual
    association between maternal education and child GCSE scores
    adjusting for the child’s polygenic score). Note that this estimate
    should be the same as a regression of Y on X adjusting for the
    polygenic score in the dataset from which the correlations were
    obtained.

  - **Genetic confounding** This is the estimate of genetic confounding.

  - **Total effect** This is the total effect which should add up to the
    observed initial association between X and Y.

## Heritability scenario

We will now implement the sensitivity analysis to examine genetic
confounding under a scenario in which polygenic scores explain
SNP-heritability in the outcome (here child GCSE scores).

We will do this by adding a “h2” option which provides the chosen
heritability estimate.

Here, h2 = 0.31, which corresponds to the SNP-heritability of Y (child
GCSE scores).

**Heritability and genetic correlation under different scenarios**

|                              | Education |
| ---------------------------- | --------- |
| Best-Fitting Polygenic score | 0.119     |
| SNP-based scenario           | **0.31**  |
| Twin scenario                | 0.63      |

``` r
gsensY(rxy=0.3975,
       rgx = 0.2909,
       rgy = 0.3462,
       n=3785,
       h2=0.31)
```

    ##                       est    se      z      pvalue ci.lower ci.upper
    ## Adjusted Bxy        0.175 0.024  7.359  1.8564e-13    0.129    0.222
    ## Genetic confounding 0.222 0.014 15.914  5.0431e-57    0.195    0.249
    ## Total effect        0.397 0.016 25.212 2.9876e-140    0.367    0.428

The results show that under a SNP heritability scenario, the effect of
maternal education on child educational achievement is attenuated
(B=0.176) relative to when controlling for observed polygenic scores
(B=0.325).

As noted in the manuscript, it is possible to fix the ratio k between
rgy and rgx when a priori knowledge is available, for example the
genetic relationship between the child and the mother.

We do this below by specifying that rgx (i.e., the correlation between
the child’s genetics with maternal education) is half of the correlation
between the child’s genetics with their own educational achievement
(i.e., the SNP heritability).

We also specify that the rgy (i.e., the correlation between the child’s
genetics with their own educational achievement) reflects SNP
heritability (by taking the square root the heritability estimate to get
the correlation).

``` r
gsensY(rxy=0.3975,
       rgx = 0.5*sqrt(.31),
       rgy = sqrt(.31),
       n = 3785,
       h2 = 0.31)
```

    ##                       est    se      z      pvalue ci.lower ci.upper
    ## Adjusted Bxy        0.263 0.016 16.941  2.2374e-64    0.232    0.293
    ## Genetic confounding 0.135 0.006 21.208 8.1051e-100    0.122    0.147
    ## Total effect        0.397 0.015 26.092 4.5077e-150    0.368    0.427

In this fixed solution, rgx and rgy and are therefore set to their value
under the heritability scenario, rather than the observed value from the
polygenic score.

# Two polygenic scores <a name="link3"></a>

When a polygenic score is available for each of X and Y, both can be
modelled using gsensXY, for example, for maternal years of education and
BMI. The gsensXY function takes a number of additional arguments:

  - **rxy** = as before, the observed phenotypic correlations between
    study variables (here: correlation between maternal years of
    education and BMI);

  - **rg1x** = the correlations between the observed phenotype x
    (maternal education) and the observed polygenic score for education
    (adjusted for sex, age, and principal components);

  - **rg2x** = the correlations between the observed phenotype x
    (maternal education) and the observed polygenic score for BMI
    (adjusted for sex, age, and principal components);

  - **rg1y** = the correlations between the observed phenotype y (BMI)
    and the observed polygenic score for education (adjusted for sex,
    age, and principal components);

  - **rg2y** = the correlations between the observed phenotype y (BMI)
    and the observed polygenic score for BMI (adjusted for sex, age, and
    principal components);

  - **rg1g2** = the correlations between the observed polygenic scores
    for education and BMI;

  - **h2.x** = the variance explained in the observed outcome x
    (maternal education) by the observed polygenic score for education
    (**h2.x = rg1x^2**);

  - **h2.y** = the variance explained in the observed outcome y (BMI) by
    the observed polygenic score for BMI (**h2.y =rg2y^2**);

  - **print = T** = enables the examination of model parameters

## Observed polygenic scores

``` r
gsensXY(rxy=-0.0894,
              rg1x=0.2909,
              rg2x=-0.0839,
              rg1y=-0.0269,
              rg2y=0.2524,
              rg1g2=-0.1856,
              n=3663,
              h2.x=0.2909^2,
              h2.y=0.2524^2,
              print=T)
```

    ## lavaan 0.6-7 ended normally after 70 iterations
    ## 
    ##   Estimator                                        GLS
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         12
    ##                                                       
    ##   Number of observations                          3663
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                                 0.000
    ##   Degrees of freedom                                 0
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Expected
    ##   Information saturated (h1) model          Structured
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   GG1 =~                                              
    ##     G1       (lg1)    1.000    0.059   16.903    0.000
    ##   GG2 =~                                              
    ##     G2       (lg2)    1.000    0.068   14.809    0.000
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   Y ~                                                 
    ##     X        (bxy)   -0.081    0.017   -4.820    0.000
    ##     GG1     (bg1y)    0.044    0.018    2.404    0.016
    ##     GG2     (bg2y)    0.254    0.004   67.390    0.000
    ##   X ~                                                 
    ##     GG1     (bg1x)    0.285    0.003   89.869    0.000
    ##     GG2     (bg2x)   -0.031    0.017   -1.805    0.071
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   GG1 ~~                                              
    ##     GG2     (bg12)   -0.186    0.022   -8.408    0.000
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     GG1               1.000                           
    ##     GG2               1.000                           
    ##    .Y         (vy)    0.930    0.023   39.974    0.000
    ##    .X         (vx)    0.914    0.023   39.407    0.000
    ##    .G1       (vg1)    0.000    0.111    0.000    1.000
    ##    .G2       (vg2)    0.000    0.129    0.000    1.000
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     conf             -0.009    0.006   -1.343    0.179
    ##     total            -0.089    0.016   -5.431    0.000
    ## 
    ## Constraints:
    ##                                                |Slack|
    ##     bg1x+bg1g2*bg2x - (sqrt(0.08462281))         0.000
    ##     bg2y+(bg2x+bg1g2*b1)*+12*1-((0.06370576))    0.000

    ##                        est    se      z     pvalue ci.lower ci.upper
    ## Adjusted Bxy        -0.081 0.017 -4.820 1.4336e-06   -0.114   -0.048
    ## Genetic confounding -0.009 0.006 -1.343    0.17943   -0.021    0.004
    ## Total effect        -0.089 0.016 -5.431 5.5933e-08   -0.122   -0.057

</br>

Similar to the one polygenic score case the model is specified with
paths between variables X and Y and polygenic scores for X (g1) and Y
(g2). In addition to the negative variance, several parameters have
unlikely values. For example, the cross path bg1y flips from a negative
to a positive value, which would imply that a higher polygenic score for
education is linked to higher BMI. In such cases, constraints can be
imposed on the the model in the following way.

``` r
gsensXY(rxy=-0.0894,
              rg1x=0.2909,
              rg2x=-0.0839,
              rg1y=-0.0269,
              rg2y=0.2524,
              rg1g2=-0.1856,
              n=3663,
              h2.x=0.2909^2,
              h2.y=0.2524^2,
              print=T, 
              constrain='lg1 < 1 \n lg2 < 1 \n  vg1 > 0  \n vg2 > 0 \n bg1y < 0 \n bg2x < 0 \n bxy < 0')
```

    ## lavaan 0.6-7 ended normally after 196 iterations
    ## 
    ##   Estimator                                        GLS
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         12
    ##   Number of inequality constraints                   7
    ##                                                       
    ##   Number of observations                          3663
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                                 6.665
    ##   Degrees of freedom                                 0
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Expected
    ##   Information saturated (h1) model          Structured
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   GG1 =~                                              
    ##     G1       (lg1)    0.998    0.011   87.425    0.000
    ##   GG2 =~                                              
    ##     G2       (lg2)    1.000                           
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   Y ~                                                 
    ##     X        (bxy)   -0.069    0.016   -4.297    0.000
    ##     GG1     (bg1y)   -0.000       NA                  
    ##     GG2     (bg2y)    0.247    0.002  141.325    0.000
    ##   X ~                                                 
    ##     GG1     (bg1x)    0.285    0.003   94.417    0.000
    ##     GG2     (bg2x)   -0.031    0.016   -1.921    0.055
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   GG1 ~~                                              
    ##     GG2     (bg12)   -0.186    0.016  -11.753    0.000
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     GG1               1.000                           
    ##     GG2               1.000                           
    ##    .Y         (vy)    0.928    0.022   42.712    0.000
    ##    .X         (vx)    0.914    0.021   42.778    0.000
    ##    .G1       (vg1)    0.000       NA                  
    ##    .G2       (vg2)    0.000                           
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     conf             -0.021    0.004   -5.328    0.000
    ##     total            -0.089    0.016   -5.458    0.000
    ## 
    ## Constraints:
    ##                                                |Slack|
    ##     bg1x+bg1g2*bg2x - (sqrt(0.08462281))         0.000
    ##     bg2y+(bg2x+bg1g2*b1)*+12*1-((0.06370576))    0.000
    ##     1 - (lg1)                                    0.002
    ##     1 - (lg2)                                    0.000
    ##     vg1 - 0                                      0.000
    ##     vg2 - 0                                      0.000
    ##     0 - (bg1y)                                   0.000
    ##     0 - (bg2x)                                   0.031
    ##     0 - (bxy)                                    0.069

    ##                        est    se      z     pvalue ci.lower ci.upper
    ## Adjusted Bxy        -0.069 0.016 -4.297 1.7339e-05   -0.100   -0.037
    ## Genetic confounding -0.021 0.004 -5.328 9.9505e-08   -0.028   -0.013
    ## Total effect        -0.089 0.016 -5.458 4.8156e-08   -0.122   -0.057

Loadings for the latent part of the model are constrained to be less
than 1 and residual variances to be positive. Cross paths and the
residual association are constrained to be negative, ie not flip sign.

Note that constraints need to be imposed cautiously on the model for two
reasons:

1)  To avoid nonsensical constraints, for example to keep the negative
    constraints on the cross paths imposed in this example in another
    situation where the polygenic scores are expected to be positively
    associated with both X and Y.

2)  Constraining parameters can decrease standard error even when the
    value of the parameter is not changed. That is, if lg1 is estimated
    to 1, constraining it to 1 will not change the estimates but will
    prevent the model from estimating a standard error for lg1 and can
    reduce standard errors of other parameters including the estimates
    of interest.

Additional remarks regarding constraints:

  - Constraining vg1 \> 0 did not prevent the warning regarding negative
    variance, which is estimated to -0.000. This can be solved by using
    vg1 \> 1e-06. However, such negative variances of -0.000 can be left
    alone as constraining them does not change estimate values but can
    reduce standard errors.

  - All model parameters should be systematically checked after fitting
    a model, and before and after constraints are imposed. In addition
    to implausible parameters, so-called ‘heywood cases’ should be
    checked. Here, as the model is standardized, no parameter should be
    above 1 or below -1. For example, a heritability parameter above 1
    would correspond to a heritability above 100%.

## Check with one polygenic score

Findings for the constrained vs the unconstrained model diverge. We can
use gsensY to and just for the outcome-related polygenic score only,
here BMI. Note that in theory, if the polygenic score for the outcome
was perfect, adjusting for that polygenic score would be sufficient as
it would capture shared genetic influences between X and Y. Findings
from adjusting for the BMI polygenic score only converge with the
constrained version of the two polygenic score approach.

``` r
gsensY(rxy=-0.0894,
             rgx = -0.0839,
             rgy = 0.2524,
             n=3663,
             h2=0.2524^2)
```

    ##                        est    se      z     pvalue ci.lower ci.upper
    ## Adjusted Bxy        -0.069 0.016 -4.233 2.3022e-05   -0.101   -0.037
    ## Genetic confounding -0.021 0.004 -5.063 4.1229e-07   -0.029   -0.013
    ## Total effect        -0.089 0.016 -5.429 5.6544e-08   -0.122   -0.057

Note here that the vg is -0.000 here too, which does not affect the
model so the warning can be safely ignored.

## Heritability scenario

``` r
gsensXY(rxy=-0.0894,
              rg1x=0.2909,
              rg2x=-0.0839,
              rg1y=-0.0269,
              rg2y=0.2524,
              rg1g2=-0.1856,
              n=3785,
              h2.x=0.25*0.31,
              h2.y=0.186,
        constrain='lg1 < 1 \n lg2 < 1 \n  vg1 > 0  \n vg2 > 0 \n bg1y < 0 \n bg2x < 0 \n bxy < 0')
```

    ## Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
    ## variances are negative

    ##                        est    se      z     pvalue ci.lower ci.upper
    ## Adjusted Bxy        -0.043 0.017 -2.464   0.013754   -0.077   -0.009
    ## Genetic confounding -0.052 0.010 -5.193 2.0673e-07   -0.072   -0.033
    ## Total effect        -0.095 0.016 -5.832  5.484e-09   -0.127   -0.063

Heritability scenarios of interest can be modelled we two polygenic
scores by replacing h2.x and h2.y by the chosen values, here
SNP-heritability estimates. Not that in h2.x=0.25x0.31, the factor 0.25
corresponds to the genetic relatedness between child and mother
(i.e. sqrt(h2.x) is computed in the model leading to the value of the
path equal to 0.5xsqrt(0.31).)
