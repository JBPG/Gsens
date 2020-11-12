
\name{gsensXY}
\alias{gsensXY}
\title{Sensivity analysis based on two polygenic scores (X and Y)}

\usage{
gsensXY(rxy,rg1x,rg2y,rg1y,rg2x,rg1g2,n,h2.x,h2.y,constrain=NULL,print=FALSE)
}


\arguments{
   \item{rxy}{the observed phenotypic correlations between study variables}
   \item{rg1x}{the correlations between the observed phenotype X and the observed polygenic score for X}
   \item{rg2x}{the correlations between the observed phenotype X and the observed polygenic score for Y}
   \item{rg1y}{the correlations between the observed phenotype Y and the observed polygenic score for X}
   \item{rg2x}{the correlations between the observed phenotype Y and the observed polygenic score for Y}
   \item{rg1g2}{the correlations between the two observed polygenic scores}
   \item{h2.x}{the variance explained in the observed outcome X by the observed polygenic score X (h2.x = rg1x^2)}
   \item{h2.y}{the variance explained in the observed outcome Y by the observed polygenic score X (h2.y = rg2x^2)}
   \item{n}{optional parameter, enables the examination of model parameters (default = FALSE)}
   \item{constrain}{agrument to contrain model parameters (default = NULL)}
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

# Model without constrains:

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


# Model with constrains:

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

}
