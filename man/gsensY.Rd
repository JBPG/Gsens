
\name{gsensY}
\alias{gsensY}
\title{Sensivity analysis based on one polygenic score for the outcome (Y)}

\usage{
gsensY(rxy,rgx,rgy,n,h2,constrain=NULL,print=FALSE)
}


\arguments{
   \item{rxy}{the observed phenotypic correlation between exposure X and outcome Y}
   \item{rgx}{the observed correlation between the polygenic score for Y and exposure X}
   \item{rgy}{the observed correlation between the polygenic score for Y and outcome Y}
   \item{n}{sample size}
   \item{h2}{the additive genetic variance explained in outcome Y under the scenario of interest}
    \item{print}{optional, enables the examination of model parameters (default = FALSE)}
   \item{constrain}{optional, agrument to constrain model parameters (default = NULL)}
}

\value{
   The function returns a dataframe with 3 estimates
   \item{Adjusted Bxy}{the standardized estimate of the relationship between X and Y, adjusted for G}
   \item{Genetic confounding}{the estimate of genetic confounding}
   \item{Total effect}{This is the total effect}

}

\references{

1. Pingault, J.-B., O’Reilly, P. F., Schoeler, T., Ploubidis, G. B., Rijsdijk, F., & Dudbridge, F. (2018). Using genetic data to strengthen causal inference in observational research. Nature Reviews Genetics, 19(9), 566–580. https://doi.org/10.1038/s41576-018-0020-3

2. Pingault, J.-B., Rijsdijk, F., Schoeler, T., Choi, S. W., Selzam, S., Kraphol, E., O’Reilly, P. F., & Dudbridge, F. Genetic sensitivity analysis: adjusting for genetic confounding in epidemiological associations. BioRxiv.
}


\examples{

gsensY(rxy=0.3975,
             rgx = 0.2909,
             rgy = 0.3462,
             n=3785,
             h2=0.3462^2)

}
