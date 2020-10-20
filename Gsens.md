Gsens
================
Jean-Baptiste Pingault & Frank Dudbridge
20 August 2020

### Overview

#### [1: Introduction](#link1)

#### [2: One polygenic score case](#link2)

#### [3: Two polygenic scores](#link3)

Introduction<a name="link1"></a>
================================

The following script provides an example on how to run Gsens, a genetically informed sensitivity analysis, in R. Please cite the following papers where the concept was first proposed and then implemented:
1. Pingault, J.-B., O’Reilly, P. F., Schoeler, T., Ploubidis, G. B., Rijsdijk, F., & Dudbridge, F. (2018). Using genetic data to strengthen causal inference in observational research. Nature Reviews Genetics, 19(9), 566–580. <https://doi.org/10.1038/s41576-018-0020-3>
2. (UPDATE) Pingault, J.-B., Rijsdijk, F., Schoeler, T., Choi, S. W., Kraphol, E., O’Reilly, P. F., & Dudbridge, F. (In Prep). Estimating the sensitivity of associations between risk factors and outcomes to shared genetic effects. BioRxiv.

The following examples are based on the correlation matrix between polygenic scores and variables, which is available as eTable 1 in the supplementary material of article 2.

To begin, please download the "gens.source.R" file from [here](https://github.com/JBPG/Gsens/blob/master/gsens.source.R) and save it into your current working directory. Then load the Gsens source code with the code below:

``` r
HOME=getwd()
setwd(HOME)
source('gsens.source.R')
```

You should now see three functions in your local environment:
- gsensY: sensivity analysis based on one polygenic score for the outcome (Y).
- gsensX: sensivity analysis based on one polygenic score for the exposure (X).
- gsensXY: sensivity analysis based on two polygenic scores for X and Y.

One polygenic score case<a name="link2"></a>
============================================

Observed scenario
-----------------

In the following example, we will test the association between based on maternal years of education (X) and child GCSE scores (Y) after controlling for the best fitting polygenic score for years of education estimated in the child.

To do so, we will use the "gsensY" function.

For this, we need to specify 5 things:

-   rxz (the correlation between the exposure \[x\] and outcome \[y\])\*
-   rgx (the correlation between the polygenic score for Y with the exposure \[x\])\*
-   rgy (the correlation between the polygenic score for Y with the outcome \[y\])\*
-   n (sample size)
-   h2 (the variance explained in the outcome, here by the observed polygenic score \[hence why h2 is rgy^2\])

\*Note that rxz, rgx, and rgy should be partial correlations adjusting for sex, age, and principal components

``` r
round(gsensY(rxy=0.3975,rgx = 0.2894,rgy = 0.3446,n=3785,h2=0.3446^2),3)
```

    ##                       est    se      z pvalue ci.lower ci.upper
    ## Adjusted Bxy        0.325 0.017 19.467      0    0.292    0.358
    ## Genetic confounding 0.073 0.003 22.840      0    0.066    0.079
    ## Total effect        0.397 0.015 26.483      0    0.368    0.427

The output provides:

-   Adjusted Bxy: This is the standardized estimate of the relationship between X and Y, adjusted for G (i.e. the residual association between maternal education and child GCSE scores adjusting for the child's polygenic score). Note that this estimate should be the same as a regression of Y on X adjusting for the polygenic score in the dataset from which the correlations were obtained.
-   Genetic confounding: This is the estimate of genetic confounding.
-   Total effect: This is the total effect which should add up to the observed initial association between X and Y.

Heritability scenario
---------------------

We will now implement the sensitivity analysis to examine genetic confounding under a scenario in which polygenic scores explain SNP-heritability in the outcome (here child GCSE scores).

We will do this by adding a "h2" option which provides the chosen heritability estimate.

Here, h2 = 0.31, which corresponds to the SNP-heritability of Y (child GCSE scores).

``` r
round(gsensY(rxy=0.3975,rgx = 0.2894,rgy = 0.3446,n=3785,h2=0.31),3)
```

    ##                       est    se      z pvalue ci.lower ci.upper
    ## Adjusted Bxy        0.176 0.024  7.357      0    0.129    0.222
    ## Genetic confounding 0.222 0.014 15.849      0    0.195    0.249
    ## Total effect        0.398 0.016 25.213      0    0.367    0.428

The results show that under a SNP heritability scenario, the effect of maternal education on child educational achievement is attenuated (B=0.176) relative to when controlling for observed polygenic scores (B=0.325).

As noted in the manuscript, it is possible to fix the ratio k between rgy and rgx when a priori knowledge is available, for example the genetic relationship between the child and the mother.

We do this below by specifying that rgx (i.e., the correlation between the child's genetics with maternal education) is half of the correlation between the child's genetics with their own educational achievement (i.e., the SNP heritability).

We also specify that the rgy (i.e., the correlation between the child's genetics with their own educational achievement) reflects SNP heritability (by taking the square root the heritability estimate to get the correlation).

``` r
 round(gsensY(rxy=0.3975,rgx = 0.5*sqrt(.31),rgy = sqrt(.31),n=3785,h2=0.31),3)
```

    ##                       est    se      z pvalue ci.lower ci.upper
    ## Adjusted Bxy        0.263 0.016 16.941      0    0.232    0.293
    ## Genetic confounding 0.135 0.006 21.208      0    0.122    0.147
    ## Total effect        0.397 0.015 26.092      0    0.368    0.427

In this fixed solution, rgx and rgy and are therefore set to their value under the heritability scenario, rather than the observed value from the polygenic score.

Two polygenic scores <a name="link3"></a>
=========================================

When polygenic scores are available for both X and Y, both can be modelled using gsensXY, for example, for maternal years of education and BMI.

Observed polygenic scores
-------------------------

``` r
round(gsensXY(rxy=-0.0894,rg1x=0.2894,rg2y=0.2522,rg1y=-0.0268,rg2x=-0.0837,rg1g2=-0.1847,n=3663,h2.x=0.2894^2,h2.y=0.2522^2,print=T),3)
```

    ## Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
    ## variances are negative

    ## lavaan 0.6-6 ended normally after 73 iterations
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
    ##     G1       (lg1)    1.000    0.059   16.823    0.000
    ##   GG2 =~                                              
    ##     G2       (lg2)    1.000    0.068   14.798    0.000
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   Y ~                                                 
    ##     X        (bxy)   -0.081    0.017   -4.817    0.000
    ##     GG1     (bg1y)    0.043    0.018    2.389    0.017
    ##     GG2     (bg2y)    0.253    0.004   67.675    0.000
    ##   X ~                                                 
    ##     GG1     (bg1x)    0.284    0.003   89.799    0.000
    ##     GG2     (bg2x)   -0.031    0.017   -1.825    0.068
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   GG1 ~~                                              
    ##     GG2     (bg12)   -0.185    0.022   -8.375    0.000
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     GG1               1.000                           
    ##     GG2               1.000                           
    ##    .Y         (vy)    0.930    0.023   39.981    0.000
    ##    .X         (vx)    0.915    0.023   39.445    0.000
    ##    .G1       (vg1)   -0.000    0.112   -0.000    1.000
    ##    .G2       (vg2)   -0.000    0.129   -0.000    1.000
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     conf             -0.009    0.006   -1.361    0.173
    ##     total            -0.089    0.016   -5.431    0.000
    ## 
    ## Constraints:
    ##                                                |Slack|
    ##     bg1x+bg1g2*bg2x - (sqrt(0.08375236))         0.000
    ##     bg2y+(bg2x+bg1g2*b1)*+12*1-((0.06360484))    0.000

    ##                        est    se      z pvalue ci.lower ci.upper
    ## Adjusted Bxy        -0.081 0.017 -4.817  0.000   -0.114   -0.048
    ## Genetic confounding -0.009 0.006 -1.361  0.173   -0.021    0.004
    ## Total effect        -0.089 0.016 -5.431  0.000   -0.122   -0.057

Similar to the one polygenic score case the model is specified with paths between variables X and Y and polygenic scores for X (g1) and Y (g2). The option print = T enables the examination of model parameters. In addition to the negative variance, several parameters have unlikely values. For example, the cross path bg1y flips from a negative to a positive value, which would imply that a higher polygenic score for education is linked to higher BMI. In such cases, constraints can be imposed on the the model in the following way.

``` r
round(gsensXY(rxy=-0.0894,rg1x=0.2894,rg2y=0.2522,rg1y=-0.0268,rg2x=-0.0837,rg1g2=-0.1847,n=3663,h2.x=0.2894^2,h2.y=0.2522^2,print=T, constrain='lg1 < 1 \n lg2 < 1 \n  vg1 > 0  \n vg2 > 0 \n bg1y < 0 \n bg2x < 0 \n bxy < 0'),3)
```

    ## Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
    ## variances are negative

    ## lavaan 0.6-6 ended normally after 191 iterations
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
    ##   Test statistic                                 6.574
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
    ##     G1       (lg1)    0.998    0.011   87.404    0.000
    ##   GG2 =~                                              
    ##     G2       (lg2)    1.000       NA                  
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   Y ~                                                 
    ##     X        (bxy)   -0.069    0.016   -4.301    0.000
    ##     GG1     (bg1y)    0.000       NA                  
    ##     GG2     (bg2y)    0.246    0.002  141.351    0.000
    ##   X ~                                                 
    ##     GG1     (bg1x)    0.284    0.003   94.300    0.000
    ##     GG2     (bg2x)   -0.031    0.016   -1.942    0.052
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   GG1 ~~                                              
    ##     GG2     (bg12)   -0.185    0.016  -11.691    0.000
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     GG1               1.000                           
    ##     GG2               1.000                           
    ##    .Y         (vy)    0.928    0.022   42.713    0.000
    ##    .X         (vx)    0.915    0.021   42.779    0.000
    ##    .G1       (vg1)   -0.000       NA                  
    ##    .G2       (vg2)   -0.000       NA                  
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     conf             -0.021    0.004   -5.314    0.000
    ##     total            -0.089    0.016   -5.458    0.000
    ## 
    ## Constraints:
    ##                                                |Slack|
    ##     bg1x+bg1g2*bg2x - (sqrt(0.08375236))         0.000
    ##     bg2y+(bg2x+bg1g2*b1)*+12*1-((0.06360484))    0.000
    ##     1 - (lg1)                                    0.002
    ##     1 - (lg2)                                    0.000
    ##     vg1 - 0                                      0.000
    ##     vg2 - 0                                      0.000
    ##     0 - (bg1y)                                   0.000
    ##     0 - (bg2x)                                   0.031
    ##     0 - (bxy)                                    0.069

    ##                        est    se      z pvalue ci.lower ci.upper
    ## Adjusted Bxy        -0.069 0.016 -4.301      0   -0.100   -0.037
    ## Genetic confounding -0.021 0.004 -5.314      0   -0.028   -0.013
    ## Total effect        -0.089 0.016 -5.458      0   -0.122   -0.057

Loadings for the latent part of the model are constrained to to be less than 1 and residual variances to be positive. Cross paths and the residual association are constrained to be negative, ie not flip sign.

Note that constraints need to be imposed cautiously on the model for two reasons:

-To avoid nonsensical constraints, for example to keep the negative constraints on the cross paths imposed in this example in another situation where the polygenic scores are expected to be positively associated with both X and Y.

-   Constraining parameters can decrease standard error even when the value of the parameter is not changed. That is, if lg1 is estimated to 1, constraining it to 1 will not change estimates but will prevent the model from estimating a standard error for lg1 and can reduce standard errors of other parameters including the estimates of interest.

Additional remarks regarding constraints:
-Constraining vg1 &gt; 0 did not prevent the warning regarding negative variance, which is estimated to -0.000. This can be solved by using vg1 &gt; 1e-06. However, such negative variances of -0.000 can be left alone as constraining them does not change estimate values but can reduce standard errors. -All model parameters should be systematically checked after fitting a model, and before and after constraints are imposed. In addition to unplausible parameters, so-called 'heywood cases' should be checked. here, as the model is standardized, no parameter should be above 1 or below -1. For example, a heritability parameter above 1 would correspond to a heritability above 100%.

Check with one polygenic score
------------------------------

Findings for the constrained vs the unconstrained model diverge. We can use gsensY to and just for the outcome-related polygenic score only, here BMI. Note that in theory, if the polygenic score for the outcome was perfect, adjusting for that polygenic score would be sufficient as it would capture shared genetic influences between X and Y. Findings from adjusting for the BMI polygenic score only converge with the constrained version of the two polygenic score appraoch.

``` r
round(gsensY(rxy=-0.0894,rgx = -0.0837,rgy = 0.2522,n=3663,h2=0.2522^2),3)
```

    ## Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
    ## variances are negative

    ##                        est    se      z pvalue ci.lower ci.upper
    ## Adjusted Bxy        -0.069 0.016 -4.238      0   -0.101   -0.037
    ## Genetic confounding -0.021 0.004 -5.052      0   -0.029   -0.013
    ## Total effect        -0.089 0.016 -5.429      0   -0.122   -0.057

Note here that the vg is -0.000 here too, which does not affect the model so the warning can be safely ignored.

Heritability scenario
---------------------

``` r
round(gsensXY(rxy=-0.0894,rg1x=0.2894,rg2y=0.2522,rg1y=-0.0268,rg2x=-0.0837,rg1g2=-0.1847,n=3785,h2.x=0.25*0.31,h2.y=0.186,
        constrain='lg1 < 1 \n lg2 < 1 \n  vg1 > 0  \n vg2 > 0 \n bg1y < 0 \n bg2x < 0 \n bxy < 0'),3)
```

    ## Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
    ## variances are negative

    ##                        est    se      z pvalue ci.lower ci.upper
    ## Adjusted Bxy        -0.043 0.017 -2.447  0.014   -0.077   -0.008
    ## Genetic confounding -0.052 0.010 -5.199  0.000   -0.072   -0.033
    ## Total effect        -0.095 0.016 -5.829  0.000   -0.127   -0.063

Heritability scenarios of interest can be modelled we two polygenic scores by replacing h2.x and h2.y by the chosen values, here SNP-heritability estimates. Not that in h2.x=0.25x0.31, the factor 0.25 corresponds to the genetic relatedness between child and mother (i.e. sqrt(h2.x) is computed in the model leading to the value of the path equal to 0.5xsqrt(0.31).)
