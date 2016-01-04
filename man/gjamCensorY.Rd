\name{gjamCensorY}
\alias{gjamCensorY}
\title{Censor gjam data}
\description{
  Returns a list with censored values, intervals, and censored response data \code{y}.
}
\usage{
  gjamCensorY(values, intervals, y, whichcol = c(1:ncol(y)))
}
\arguments{
  \item{values}{Values in \code{y} that are censored, specified by \code{intervals}}

  \item{intervals}{Matrix having two rows and one column for each value 
  in \code{values}. The first row holds lower bounds.  The second row holds
  upper bounds.  See \bold{Examples}.}

  \item{y}{Response data matrix, \code{n} rows by \code{S} columns.  All values within \code{intervals} will be replaced with \code{values}}

  \item{whichcol}{Columns in \code{y} that are censored (often not all responses are
  censored)}
}
\value{
  Returns a list containing two elements.

  \item{\code{y}}{n by S matrix updated with censored {values} substituted for
  those falling within \code{intervals}.}

  \item{\code{censor}}{list containing \code{$columns} that are censored and
  \code{$partition}, a matrix with 3 rows used in \code{\link{gjamGibbs}} 
  and \code{\link{gjamPlot}}, one column per censor interval.  Rows 
  are values, followed by lower and upper bounds.}
}
\details{
  Any values in \code{y} that fall within censored \code{intervals} are replaced with
  censored \code{values}.
}
\examples{
## trait data for USFS For Inventory $ Analysis program
data(forestTraits)
y <- forestTraits$y

# continuous traits censored to allow negative to positive range
wc       <- 1:6
tmp      <- gjamCensorY(values=-1000,intervals=cbind( c(-Inf,-1000) ),
                         y = y,whichcol=wc)
censor <- list('CA' = tmp$censor)

# continuous traits on (0, 1)
wc     <- 13:14
tmp    <- gjamCensorY(values=c(0,1),intervals=cbind( c(-Inf,0),c(1,Inf) ),
                       y = y,whichcol=wc)
censor <- append(censor,list('CA' = tmp$censor))

# an example with a few iterations (ng = 200)
modelList <- list(ng=20,burnin=10, typeNames = forestTraits$typeNames,
                  xfactors='soilFactor',censor=censor)
outTrait  <- gjamGibbs(~temp+moisture*(deficit+soilFactor),
                        xdata = forestTraits$xdata, y = y, 
                        modelList = modelList)
                        
# repeat with ng = 3000, burnin = 500, then plot data here:
plotPars  <- list(width=3,height=2)                  
fit       <- gjamPlot(output = outTrait, plotPars)
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
  \code{\link{gjamGibbs}} analyzes data
  
  A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

The manual and data are available at 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

