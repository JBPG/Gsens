
\name{gsensXY}
\alias{gsensXY}
\title{Sensivity analysis based on two polygenic scores (X and Y)}

\usage{
gsensXY(rxy,rg1x,rg2y,rg1y,rg2x,rg1g2,n,h2.x,h2.y,constrain=NULL,print=FALSE)
}


\arguments{
   \item{rxy}{the observed phenotypic correlation between exposure X and outcome Y}
   \item{rg1x}{the correlation between the exposure X and the observed polygenic score for X}
   \item{rg2y}{the correlation between the outcome Y and the observed polygenic score for Y}
   \item{rg1y}{the correlation between the outcome Y and the observed polygenic score for X}
   \item{rg2x}{the correlation between the exposure X and the observed polygenic score for Y}
   \item{rg1g2}{the correlation between the two observed polygenic scores}
   \item{n}{sample size}
   \item{h2.x}{the additive genetic variance explained in exposure X under the scenario of interest}
   \item{h2.y}{the additive genetic variance explained in outcome Y under the scenario of interest}
   \item{print}{optional, enables the examination of model parameters (default = FALSE)}
   \item{constrain}{optional, agrument to constrain model parameters (default = NULL)}
}


\value{
   The function returns a dataframe with 3 estimates
   \item{Adjusted Bxy}{the standardized estimate of the relationship between X and Y, adjusted for G1 and G2}
   \item{Genetic confounding}{the estimate of genetic confounding}
   \item{Total effect}{This is the total effect}
}

\references{

1. Pingault, J.-B., O’Reilly, P. F., Schoeler, T., Ploubidis, G. B., Rijsdijk, F., & Dudbridge, F. (2018). Using genetic data to strengthen causal inference in observational research. Nature Reviews Genetics, 19(9), 566–580. https://doi.org/10.1038/s41576-018-0020-3

2. Pingault, J.-B., Rijsdijk, F., Schoeler, T., Choi, S. W., Selzam, S., Kraphol, E., O’Reilly, P. F., & Dudbridge, F. Genetic sensitivity analysis: adjusting for genetic confounding in epidemiological associations. BioRxiv.
}


\examples{
}
