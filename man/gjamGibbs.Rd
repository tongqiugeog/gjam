\name{gjamGibbs}
\alias{gjamGibbs}
\title{Gibbs sampler for gjam data}
\description{
  Analyzes joint attribute data (e.g., species abundance) with Gibbs sampling.  
  Can be output from \code{\link{gjamSimData}}.  For a data example see
  \code{\link{gjamCensorY}}.  Returns a list of objects from 
  Gibbs sampling that can be plotted by \code{\link{gjamPlot}}.
}
\usage{
  gjamGibbs(formula, xdata, y, modelList)
}
\arguments{
  \item{formula}{R formula for model, e.g., \code{~ x1 + x2}.}

  \item{xdata}{\code{data.frame} for predictors to be included in model must
  include column names that match variable names in \code{formula}.}

  \item{y}{\code{n} by \code{S} response matrix.  Column names are 
  unique attribute labels, e.g., species names.  All columns will be 
  included in the analysis.}

  \item{modelList}{List specifying prior distributions, holdouts, and 
  traits.  Must include \code{ng} (number of Gibbs steps), \code{burnin}, 
  and \code{typeNames}. Can include number of holdouts for out-of-sample 
  prediction, \code{holdoutN}.  See \bold{Details}.}
}
\details{
  Constructs a design matrix \code{x} based on columns in \code{xdata} and
  \code{formula}.  Note that \code{formula} does not contain \code{y}, which
  is passed as a \code{n} by \code{S} matrix.  Both predictors \code{x} and 
  responses \code{y} can include missing values as \code{NA}.  \code{modelList} 
  has these defaults or options:
  
  \code{ng} = 3000, number of Gibbs steps.  
    
  \code{burnin} = 1000, no. initial steps, must be less than \code{ng}.  
    
  \code{typeNames} can be \code{'PA'} (presenceAbsence), \code{'CA'}  
  (continuous), \code{'DA'} (discrete),  \code{'FC'} (fractional composition),
  \code{'CC'} (count composition), \code{'OC'} (ordinal).  \code{typeNames} can 
  be a single value to be applied to all columns in \code{y} or there can be one 
  value for each column.
    
  \code{holdoutN} = 0, number of observations to hold out for out-of-sample
    prediction.  
    
  \code{holdoutIndex} = numeric(0), vector of plot numbers to holdout for
    out-of-sample prediction.  
    
  \code{censor} = NULL, list specifying columns, values, and intervals for
    censoring, see \code{\link{gjamCensorY}}.
    
  \code{effort} = NULL, list containing \code{'columns'}, a vector of length
    \code{<= S} giving the names of columns in in \code{y}, and \code{'values'},
     a length-\code{n} vector of effort (see \bold{Examples}).  \code{effort} is plot
     area for discrete count data \code{'DA'}.
     
  \code{xfactors} = NULL, character vector indicating variables in \code{formula}
  to be treated as factors.
  
  A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

The manual and data are available at 'http://sites.nicholas.duke.edu/clarklab/code/'.

}
\value{
  \item{missingIndex}{index for missing values in design \code{x}.}

  \item{xMissMu,xMissSd}{ posterior mean and se for predicted \code{x}.}

  \item{yMissMu,yMissSd}{ posterior mean and se for predicted \code{y}.}

  \item{chains}{a list of MCMC chains, including \code{rgibbs} 
  (correlation: \code{ng} by \code{S^2}), \code{sensGibbs} 
  (\code{ng} by \code{q} sensitivity), \code{sgibbs} (\code{ng} by 
  \code{S^2} covariance), \code{bgibbs}(\code{ng} by \code{q*S} coefficients)}

  \item{richness}{ posterior estimates of no. of species present, one row per plot}

  \item{modelSummary}{list containing the following objects:}

    \item{designTable}{Variance inflation factors and correlations for 
    predictors}

    \item{classBySpec}{\code{S} by \code{K} summary of \code{y}, 
    where \code{K} is the number of of classes.}

    \item{yMu,ySd}{ predictions of \code{y} with standard errors}

    \item{betaMu,betaSe}{posterior means and standard errors for coefficients.}
    
    \item{sigMu,sigSe}{ posterior means and standard errors for covariances.}
    
    \item{corMu,corSe}{ posterior means and standard errors for correlations.}
    
    \item{cutMu,cutSe}{ posterior means and standard errors for ordinal 
    partition.}

    \item{xpredMu,xpredSd}{ inverse predictions of \code{x}.}

    \item{wMu,wSd}{ posterior estimates of latent states.}
}
\examples{
## combinations of scales
typeNames <- c('OC','OC','OC','CC','CC','CC',
               'CC','CC','CA','CA','PA','PA')         
simMIX    <- gjamSimData(n=1000,S=length(typeNames),q=3,typeNames=typeNames)
modelList <- list(ng = 50, burnin = 5, typeNames = simMIX$typeNames)
outMIX    <- gjamGibbs(simMIX$formula, simMIX$xdata, simMIX$y, modelList)

# repeat with ng = 2000, burnin = 500, then plot data:
plotPars  <- list(trueValues = simMIX$trueValues,width=3,height=2)
fit       <- gjamPlot(output = outMIX, plotPars)

## discrete abundance with heterogeneous effort 
S         <- 5                             
n         <- 1000
effort    <- list( columns = 1:S, values = round(runif(n,.5,5),1) )
simDA     <- gjamSimData(n,S,q=5,typeNames='DA',effort=effort)
modelList <- list(ng=50, burnin=5, typeNames = simDA$typeNames,
                  effort = effort)
outDA     <- gjamGibbs(simDA$formula,simDA$xdata, simDA$y, modelList)

# repeat with ng = 2000, burnin = 500, then plot data:
plotPars  <- list(trueValues = simDA$trueValues,width=3,height=2)
fit       <- gjamPlot(output = outDA, plotPars)
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. (in review) 
Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.
}
\seealso{
\code{\link{gjamSimData}} simulates data

A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

The manual and data are available at 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

