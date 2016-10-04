\name{gjam-package}
\alias{gjam-package}
\alias{gjam}
\alias{GJAM}
\docType{package}
\concept{multivariate regression}
\concept{hierarchical Bayes}
\concept{Gibbs sampling}
\concept{dimension reduction}
\concept{species distribution model}
\concept{ecological traits}
\concept{prediction}
\concept{ordination}
\title{
  Generalized Joint Attribute Modeling
}
\description{
  Inference and prediction for jointly distributed responses that 
  are combinations of continous and discrete data.  Functions begin with \code{'gjam'} to avoid conflicts with other packages.
}
\details{
  \tabular{ll}{
  Package: \tab gjam\cr
  Type: \tab Package\cr
  Version: \tab 2.0\cr
  Date: \tab 2015-12-21\cr
  License: \tab GPL (>= 2)\cr
  URL: \tab http://sites.nicholas.duke.edu/clarklab/code/\cr
  }
  
The generalized joint attribute model (gjam) analyzes multivariate data that are combinations of presence-absence, ordinal, continuous, discrete, composition, zero-inflated, and censored.  It does so as a joint distribution over response variables.  gjam provides inference on sensitivity to input variables, correlations between responses on the data scale, model selection, and prediction.  

Importantly, analysis is done on the observation scale.  That is, coefficients and covariances are interpreted on the same scale as the data.  Contrast this approach with standard Generalized Linear Models, where coefficients and covariances are difficult to interpret and cannot be compared across responses that are modeled on different scales.

gjam was motivated by species distribution and abundance data in ecology, but can provide an attractive alternative to traditional methods wherever observations are multivariate and combine multiple scales and mixtures of continuous and discrete data.  

gjam can be used to model ecological trait data, where species traits are translated to locations as community-weighted means and modes.

Posterior simulation is done by Gibbs sampling.  Analysis is done by these functions:


\code{\link{gjamSimData}} simulates data for analysis by \code{gjamGibbs}.

\code{\link{gjamCensorY}} defines censored values and intervals.

\code{\link{gjamTrimY}} trims the response matrix and aggregates rare types.

\code{\link{gjamGibbs}} fits model with Gibbs sampling.

\code{\link{gjamPlot}} plots output from \code{gjamGibbs}.

\code{\link{gjamSpec2Trait}} constructs plot by trait matrix.

\code{\link{gjamPredict}} does conditional prediction.

\code{\link{gjamOrdination}} ordinates the response matrix.

\code{\link{gjamDeZero}} de-zeros response matrix for storage.

\code{\link{gjamReZero}} recovers response matrix from de-zeroed format.

\code{\link{gjamIIE}} evaluates indirect effects and interactions.

\code{\link{gjamIIEplot}} plots indirect effects and interactions.

\code{\link{gjamSpec2Trait}} generates trait values.

\code{\link{gjamPoints2Grid}} aggregates incidence data to counts on a lattice.

}


\author{
  Author: James S Clark, \email{jimclark@duke.edu}, Daniel Taylor-Rodriquez
}
\references{
  Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2016.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs, in press.

  Clark, J.S. 2016. Why species tell more about traits than traits tell us about species: Predictive models.  Ecology 97, 1979-1993.
  
  Taylor-Rodriguez, D., K. Kaufeld, E. M. Schliep, J. S. Clark, and A. E. Gelfand. 2016.  Joint species distribution modeling: dimension eduction using Dirichlet processes.  Bayesian Analysis, in press.
}
\keyword{ package }
\seealso{
  \code{\link{gjamSimData}}, 
  \code{\link{gjamCensorY}}, 
  \code{\link{gjamTrimY}}, 
  \code{\link{gjamGibbs}}, 
  \code{\link{gjamPredict}}, 
  \code{\link{gjamSpec2Trait}}, 
  \code{\link{gjamPlot}},
  \code{\link{gjamDeZero}},
  \code{\link{gjamReZero}}
  
  A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

The manual and data are available at 'http://sites.nicholas.duke.edu/clarklab/code/'.
}
