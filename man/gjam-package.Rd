\name{gjam-package}
\alias{gjam-package}
\alias{gjam}
\alias{GJAM}
\docType{package}
\concept{multivariate regression}
\concept{species distribution model}
\concept{ecological traits}
\title{
  Generalized Joint Attribute Modeling
}
\description{
  Inference and prediction for jointly distributed responses that 
  are combinations of continous and discrete data.
}
\details{
  \tabular{ll}{
  Package: \tab gjam\cr
  Type: \tab Package\cr
  Version: \tab 1.0\cr
  Date: \tab 2015-12-21\cr
  License: \tab GPL (>= 2)\cr
  URL: \tab http://sites.nicholas.duke.edu/clarklab/code/\cr
  }
  
The generalized joint distribution model (gjam) analyzes multivariate data that are combinations of presence-absence, ordinal, continuous, discrete, composition, zero-inflated, and censored.  It does so as a joint distribution over response variables.  gjam provides inference on sensitivity to input variables, correlations between responses on the data scale, model selection, and prediction.  

Importantly, analysis is done on the observation scale.  That is, coefficients and covariances are interpreted on the same scale as the data.  Contrast this approach with standard Generalized Linear Models, where coefficients and covariances are difficult to interpret and cannot be compared across responses that are modeled on different scales.

gjam was motivated by species distribution and abundance data in ecology, but can provide an attractive alternative to traditional methods wherever observations are multivariate and combine multiple scales and mixtures of continuous and discrete data.  

Posterior simulation is done by Gibbs sampling.  Analysis is done by four functions:

\code{\link{gjamSimData}} simulates data for analysis by \code{gjam_gibbs}.

\code{\link{gjamCensorY}} defines censored values and intervals.

\code{\link{gjamGibbs}} analyzes data with Gibbs sampling.

\code{\link{gjamPlot}} plots output from \code{gjam_gibbs}.

Several simulated data sets illustrate examples:

\code{\link{forestTraits}} contains ecological trait data from North American Forest Inventory data

}


\author{
  Author: James S Clark, \email{jimclark@duke.edu}
}
\references{
  Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang (in review).
  Generalized joint attribute modeling for biodiversity analysis, Median-zero, 
  multivariate, multifarious data.

  Clark, J.S. (in review). Why species tell more about traits than vice versa: 
Predictive analysis to connect species, traits, and environment.
}
\keyword{ package }
\seealso{
  \code{\link{gjamSimData}}, 
  \code{\link{gjamCensorY}}, 
  \code{\link{gjamGibbs}}, 
  \code{\link{gjamPlot}}
  A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

The manual and data are available at 'http://sites.nicholas.duke.edu/clarklab/code/'.
}
