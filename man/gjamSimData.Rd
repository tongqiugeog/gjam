\name{gjamSimData}
\alias{gjamSimData}
\title{Simulated data for gjam analysis}
\description{
  Simulates data for analysis by \code{\link{gjamGibbs}}.
}
\usage{
  gjamSimData(n = 1000, S = 10, q = 5, nmiss = 0, typeNames, effort = NULL)
}
\arguments{
  \item{n}{Sample size}

  \item{S}{Number of response variables (columns) in \code{y} should be less than
   \code{n}}

  \item{q}{Number of predictors (columns) in design matrix \code{x}
  should be \code{<< n}}

  \item{nmiss}{Number of missing values to in \code{x} should be less than \code{n}}

  \item{typeNames}{Character vector of data types, see \bold{Details}}

  \item{effort}{List containing '\code{columns}' specifying columns to which
  effort applies, and '\code{values}', a length-\code{n} vector of effort per 
  observation.}
}
\value{

  \item{formula}{R formula for model, e.g., \code{~ x1 + x2}}
  
  \item{xdata}{\code{data.frame} includes columns for predictors in
  design matrix}
  
  \item{y}{\code{n} by \code{S} response matrix}
  
    \item{w}{\code{n} by \code{S} latent states}
  
  \item{typeNames}{vector of data types corresponding to columns in \code{y}, see \bold{Details}} 
  
  \item{trueValues}{list containing true parameter values \code{beta} (regression
  coefficients), \code{sigma} (covariance matrix), \code{corSpec} (correlation 
  matrix corresponding to \code{sigma}), and \code{cuts} (partition matrix for 
  ordinal data).}
}
\details{
  Generates simulated data from parameters for analysis by \code{\link{gjamGibbs}}.  Because both
  parameters and data are stochastic, not all simulations will give good results.  Returns a list of objects for analysis by \code{\link{gjamGibbs}}.  
  \code{typeNames} can be '\code{PA}' (presenceAbsence), '\code{CA}'  
  (continuous), '\code{DA}' (discrete),  '\code{FC}' (fractional composition),
  '\code{CC}' (count composition), '\code{OC}' (ordinal counts).  If length-one, 
  the value in \code{typeNames} applies to all columns in \code{y}.  If not,
  \code{typeNames} is length-\code{S}, identifying each response by column 
  in \code{y}.
  
  A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

The manual and data are available at 'http://sites.nicholas.duke.edu/clarklab/code/'.

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
}
\examples{
## ordinal data
simOC <- gjamSimData(n=2000,S=10,q=3,typeNames='OC')

## continuous data censored at zero
simCA <- gjamSimData(n=2000,S=15,q=3,typeNames='CA')

## continuous and discrete data
typeNames <- c(rep('DA',5),rep('CA',4))
simDACA   <- gjamSimData(n=2000,S=length(typeNames),q=4,typeNames=typeNames)
                             
## composition count data
simCC <- gjamSimData(n=2000,S=8,q=5,typeNames='CC')

## discrete abundance, heterogeneous effort 
S      <- 5
n      <- 1500
effort <- list( columns = 1:S, values = round(runif(n,.5,5),1) )
simDA  <- gjamSimData(n,S,q=5,typeNames='DA',effort=effort)

## combinations of scales
typeNames <- c('OC','OC','OC','CC','CC','CC','CC','CC','CA','CA','PA','PA')
simMIX    <- gjamSimData(n=2000,S=length(typeNames),q=3,typeNames=typeNames)
}
