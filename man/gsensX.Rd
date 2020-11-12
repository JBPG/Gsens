
\name{gsensY}
\alias{gsensY}
\title{Sensivity analysis based on one polygenic score for the exposure (X)}

\usage{
gsensX(rxy, rgx,rgy,n,h2,constrain=NULL,print=FALSE)
}
gsensY

\arguments{
   \item{rxy}{the observed phenotypic correlations between study variables}
   \item{rgx}{the observed correlations between phenotype X}
   \item{rgy}{the observed correlations between phenotype Y}
   \item{n}{sample size}
   \item{h2}{the variance explained in the exposure X}
}

\value{
   The function returns a dataframe with 3 estimates
   \item{Adjusted Bxy}{the standardized estimate of the relationship between X and Y, adjusted for G}
   \item{Genetic confounding}{the estimate of genetic confounding}
   \item{Total effect}{This is the total effect which should add up to the observed initial association between X and Y}
}

\references{

1. Pingault, J.-B., O’Reilly, P. F., Schoeler, T., Ploubidis, G. B., Rijsdijk, F., & Dudbridge, F. (2018). Using genetic data to strengthen causal inference in observational research. Nature Reviews Genetics, 19(9), 566–580. https://doi.org/10.1038/s41576-018-0020-3

2. Pingault, J.-B., Rijsdijk, F., Schoeler, T., Choi, S. W., Kraphol, E., O’Reilly, P. F., & Dudbridge, F. (In Prep). Estimating the sensitivity of associations between risk factors and outcomes to shared genetic effects. BioRxiv.
}


\examples{

}
