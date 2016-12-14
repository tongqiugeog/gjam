\name{gjam}
\alias{gjam}
\alias{print.gjam}
\alias{summary.gjam}
\title{Gibbs sampler for gjam data}
\description{
  Analyzes joint attribute data (e.g., species abundance) with Gibbs sampling.  
  Input can be output from \code{\link{gjamSimData}}.  Returns a list of objects from 
  Gibbs sampling that can be plotted by \code{\link{gjamPlot}}.
}
\usage{
  gjam(formula, xdata, ydata, modelList)
  
  \method{print}{gjam}(x, ...)
  
  \method{summary}{gjam}(object, ...)
}
\arguments{
  \item{formula}{R formula for model, e.g., \code{~ x1 + x2}.}

  \item{xdata}{\code{data.frame} containing predictors in \code{formula}.  If not found in \code{xdata} variables are taken from \code{environment(formula)}.}

  \item{ydata}{\code{n} by \code{S} response \code{matrix} or \code{data.frame}.  Column names are unique labels, e.g., species names.  All columns will be included in the analysis.}

  \item{modelList}{\code{list} specifying inputs, including \code{ng} (number of Gibbs steps), \code{burnin}, and \code{typeNames}. Can include number of holdouts for out-of-sample prediction, \code{holdoutN}.  See \bold{Details}.}
  
  \item{x}{object of \code{\link{gjam}}.}
  
  \item{object}{currently, also an object of \code{\link{gjam}}.}
  
  \item{...}{further arguments not used here.}
}

\details{


Note that \code{formula} begins with \code{~}, not \code{y ~}.  The response matrix is passed in the form of a \code{n} by \code{S} matrix or \code{data.frame} \code{ydata}, passed as a separate argument.  
  
Both predictors in \code{xdata} and responses \code{ydata} can include missing values as \code{NA}. Factors in \code{xdata} should be declared using \code{factor}.   For computational stability variables that are not factors are standardized by mean and variance, then transformed back to original scales in \code{output}.  To retain a variable in its original scale during computation include it in the character string \code{notStandard} as part of the \code{list modelList}. (example shown in the \code{vignette} on traits).  
  
  
\code{modelList} has these defaults and provides these options:
  
  \code{ng = 2000}, number of Gibbs steps.  
    
  \code{burnin = 500}, no. initial steps, must be less than \code{ng}.  
    
  \code{typeNames} can be \code{'PA'} (presenceAbsence), \code{'CON'} (continuous on \code{(-Inf, Inf)}), \code{'CA'} (continuous abundance, zero censoring), \code{'DA'} (discrete abundance),  \code{'FC'} (fractional composition),
  \code{'CC'} (count composition), \code{'OC'} (ordinal counts), \code{'CAT'} (categorical classes).  \code{typeNames} can be a single value that applies to all columns in \code{ydata}, or there can be one value for each column.
    
  \code{holdoutN = 0}, number of observations to hold out for out-of-sample
    prediction.  
    
  \code{holdoutIndex = numeric(0)}, \code{numeric vector} of observations (row numbers) to holdout for out-of-sample prediction.  
    
  \code{censor = NULL}, \code{list} specifying columns, values, and intervals for
    censoring, see \code{gjamCensorY}.
    
  \code{effort = NULL}, \code{list} containing \code{'columns'}, a vector of length \code{<= S} giving the names of columns in in \code{y}, and \code{'values'}, a length-\code{n} vector of effort or a \code{n} by \code{S} matrix (see Examples).  \code{effort} can be plot area, search time, etc. for discrete count data \code{'DA'}.
     
  \code{notStandard = NULL}, \code{character vector} of columns in \code{xdata} that should not be standardized.
  
  \code{reductList = list(N = 20, r = 3)}, \code{list} of dimension reduction parameters.  See \code{vignette} on Dimension Reduction.
  
  \code{FCgroups, CCgroups}, are length-\code{S vectors} assigned to columns in \code{ydata} indicating composition \code{'FC'} or \code{'CC'} group membership.  For example, if there are two 'CA' columns in \code{ydata} followed by two groups of fractional composition data, each having three columns, then \code{typeNames = c('CA','CA','FC','FC','FC','FC','FC','FC')} and \code{FCgroups = c(0,0,1,1,1,2,2,2)}.  note: \code{gjamSimData} is not currently set up to simulate multiple composition groups, but \code{gjamGibbs} can model it.
  
  \code{PREDICTX = T} executes inverse prediction of \code{x}.  Speed-up by setting \code{PREDICTX = F}.
  
    \code{ematAlpha = .5} is the probability assigned for conditional and marginal independence in the \code{ematrix}.
  
  \code{traitList = list(plotByTrait, traitTypes, specByTrait)}, list of trait objects.  See vignette on Trait analysis.
  
  
More detailed vignettes can be obtained with:

\code{browseVignettes('gjam')}
}

\value{

Returns an object of \code{\link{class}} "gjam", which is a list containing the following components:

  \item{breakMat}{partition matrix.}
  
  \item{censor}{\code{list}, see \code{gjamCensor}. }
  
  \item{chains}{\code{list} of MCMC chains, all having \code{ng} rows; includes coefficients \code{bgibbs}(\code{Q*S} columns), sensitivity \code{fgibbs} (\code{Q1} columns), and \code{fbgibbs} (\code{Q1} columns, where \code{Q1 = Q - 1}, unless there are multilevel factors); covariance \code{sgibbs} has \code{S*(S + 1)/2} columns (\code{REDUCT == F}) or \code{N*r} columns (\code{REDUCT == T}).}
  
  \item{effort}{a \code{list} of effort columns and values, see Details.}
    
  \item{holdoutIndex}{index for holdout observations.}

  \item{missingIndex}{index for missing values in design \code{x}.}
  
  \item{modelList}{\code{list} of inputs.}
  
  \item{otherpar}{additional parameters used in dimension reduction.}
  
  \item{REDUCT}{\code{logical}, indicating if dimension reduction; see Details.}
  
  \item{standX}{mean and variance of standardized values in \code{x}.}
  
  \item{traitList}{inputs for trait analysis, see trait \code{vignette}.}
  
  \item{TRAITS}{\code{logical} indicating if trait prediction, see trait \code{vignette}.}
  
  \item{x}{design \code{matrix}.}
  
  \item{xdata}{input \code{data.frame}.}
  
  \item{xMissMu,xMissSd}{ posterior mean and se for missing \code{x}.}

  \item{yMissMu,yMissSd}{ posterior mean and se for missing \code{y}.}
  
  

 There are two lists containing output. \bold{\code{modelSummary}} is a \code{list} containing the following objects:

    \item{classBySpec}{\code{S} by \code{K} summary table of \code{ydata}, 
    where \code{K} is the number of classes.}
    
    \item{contrasts}{ if there are factors.}
    
    \item{DIC}{ deviance information criterion.}
    
    \item{ematAlpha}{ probability level for assigning independence in \code{ematrix}.}
    
    \item{factorList, isFactor}{ identify factors in the model in two formats.}

    \item{richness}{ obs and posterior estimates of no. of species present, one row per observation.}
    
    \item{RMSPE}{ root mean squared prediction error, by column.}
    
    \item{whConZero,whichZero}{ elements of \code{ematrix} that are zero, conditionally and marginally, at probability level \code{ematAlpha}}
  
    \item{xpredMu,xpredSd}{ inverse predictions of \code{x}.}

    \item{wMu,wSd}{ posterior estimates of latent states.}
    
    \item{ypredMu,ypredSd}{ predictions of \code{y} with standard errors}

    
 \bold{\code{parameterTables}} is a \code{list} containing the following objects:
  
    \item{betaMu,betaSe}{posterior means and standard errors for coefficients.}
    
    \item{betaTraitMu,betaTraitSe}{posterior means and standard errors for trait coefficients.}
    
    \item{corMu,corSe}{ posterior means and standard errors for correlations.}
    
    \item{cutMu,cutSe}{ posterior means and standard errors for ordinal 
    partition.}
    
    \item{ematrix}{ posterior mean E matrix on correlation scale.}
    
     \item{fBetaMu,fBetaSd}{ posterior means and standard errors for sensitivity coefficients, dimensionless and (if factors included) with contrasts.}
  
    \item{fmatrix}{ posterior sensitivity matrix F.}
    
    \item{fMu,fSe}{ posterior means and standard errors for sensitivity coefficients.}
    
    \item{sigMu,sigSe}{ posterior means and standard errors for covariances.}
    
   \item{sigmaTraitMu,sigmaTraitSe}{ posterior means and standard errors for trait covariances.}
   
      \item{tMu,tMuOrd,tSd}{ posterior means and standard errors for trait values, including ordinal traits.}
  
}
\examples{
## combinations of scales
types  <- c('DA','DA','OC','OC','OC','OC','CC','CC','CC','CC','CC','CA','CA','PA','PA')         
f      <- gjamSimData(S = length(types), typeNames = types)
ml     <- list(ng = 50, burnin = 5, typeNames = f$typeNames)
out    <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)

# repeat with ng = 5000, burnin = 500, then plot data:
pl  <- list(trueValues = f$trueValues)
gjamPlot(out, plotPars = pl)

## discrete abundance with heterogeneous effort 
S   <- 5                             
n   <- 1000
eff <- list( columns = 1:S, values = round(runif(n,.5,5),1) )
f   <- gjamSimData(n, S, typeNames='DA', effort=eff)
ml  <- list(ng = 50, burnin = 5, typeNames = f$typeNames, effort = eff)
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
summary(out)

# repeat with ng = 2000, burnin = 500, then plot data:
pl  <- list(trueValues = f$trueValues)
gjamPlot(out, plotPars = pl)
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2016.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs, in press.
}
\seealso{
\code{\link{gjamSimData}} simulates data

A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

