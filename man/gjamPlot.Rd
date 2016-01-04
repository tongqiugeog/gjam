\name{gjamPlot}
\alias{gjamPlot}
\title{Plot gjam analysis}
\description{
  Constructs plots of posterior distributions, predictive distributions, 
  and additional analysis from output of \code{\link{gjamGibbs}}.
}
\usage{
  gjamPlot(output, plotPars)
}
\arguments{
  \item{output}{list created by \code{\link{gjamGibbs}}}
  
  \item{plotPars}{list having default values in Details}
}
\value{
  Summary tables of parameter estimates are:
  
  \item{\code{betaEstimates}}{Posterior summary of beta coefficients.}

  \item{\code{corEstimates}}{Posterior summary of correlation matrix.}
}

\details{

  \code{plotPars} is a list that can contain the following, with default values:
  
  \tabular{ll}{
  
    \code{sigOnly} = F \tab plot distributions only for parameters having 
    credible intervals that do not contain zero. \cr
  
    \code{PLOTY} = T \tab plot predicted \code{y}.\cr
  
    \code{PLOTX} = T \tab plot inverse predicted \code{x}.\cr
  
    \code{ncluster} \tab number of clusters to highlight in 
    cluster diagrams, default based on \code{S}.\cr
  
    \code{corLines} = T \tab draw grid lines on grid plots of R and E.\cr
  
    \code{cex} = 1 \tab text size for grid plots, see \code{\link{par}}.\cr
  
    \code{betaGrid} = T \tab draw grid of beta coefficients.\cr
  
    \code{plotAllY} = F \tab an individual plot for each column in \code{y}.\cr
    
    \code{SMALLPLOTS} = T \tab avoids plot margin error on some devices, usually better appearance if \code{FALSE}.\cr
    
    \code{CLUSTERPLOTS} = F \tab cluster and grid plots derived from parameters.\cr
  
    \code{width,height} = 9 \tab can be small values, in inches, to avoid 
    plot margin error on some devices.\cr
  }
}
\examples{
## ordinal data
simOC     <- gjamSimData(n=1000,S=10,q=3,typeNames='OC') 
modelList <- list(ng = 100, burnin = 10, typeNames = simOC$typeNames,
                  holdoutN = 10)
outOC     <- gjamGibbs(simOC$formula, simOC$xdata, simOC$y, modelList)

# repeat with ng = 2000, burnin = 500, then plot data here:
plotPars  <- list(trueValues = simOC$trueValues,width=3,height=2)
fit       <- gjamPlot(output = outOC, plotPars)
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. (in review) 
Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.
}
\seealso{
\code{\link{gjamGibbs}}
A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

The manual and data are available at 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

