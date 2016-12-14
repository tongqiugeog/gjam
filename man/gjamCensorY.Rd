\name{gjamCensorY}
\alias{gjamCensorY}
\title{Censor gjam response data}
\description{
  Returns a \code{list} with censored values, intervals, and censored response \code{matrix y}.
}
\usage{
  gjamCensorY(values, intervals, y, whichcol = c(1:ncol(y)))
}
\arguments{
  \item{values}{Values in \code{y} that are censored, specified by \code{intervals}}

  \item{intervals}{\code{matrix} having two rows and one column for each value 
  in \code{values}. The first row holds lower bounds.  The second row holds
  upper bounds.  See \bold{Examples}.}

  \item{y}{Response \code{matrix}, \code{n} rows by \code{S} columns.  All values within \code{intervals} will be replaced with \code{values}}

  \item{whichcol}{Columns in \code{y} that are censored (often not all responses are
  censored)}
}
\value{
  Returns a list containing two elements.

  \item{\code{y}}{n by S matrix updated with censored {values} substituted for
  those falling within \code{intervals}.}

  \item{\code{censor}}{list containing \code{$columns} that are censored and
  \code{$partition}, a matrix with 3 rows used in \code{\link{gjam}} 
  and \code{\link{gjamPlot}}, one column per censor interval.  Rows 
  are values, followed by lower and upper bounds.}
}
\details{
  Any values in \code{y} that fall within censored \code{intervals} are replaced with censored \code{values}.  The example below simulates data collected on an 'octave scale': 0, 1, 2, 4, 8, ..., an approach to accelerate data collection with approximate bins.
}
\examples{
## continuous data censored at zero
values <- c(0, 2^c(0:4), Inf)                                # data in octaves
intvls <- rbind( c(-Inf, values[-length(values)]), values )  # intervals

f   <- gjamSimData(n = 2000, S = 15, Q = 3, typeNames='CA')
y   <- f$y
cc  <- c(3:6)                                                # censored columns
tmp <- gjamCensorY(values = values, intervals = intvls,
                         y = y, whichcol = cc)
y[,cc] <- tmp$y                                              # replace columns
censor <- list('CA' = tmp$censor)

ml     <- list(ng = 100, burnin = 20, censor = censor, typeNames = f$typeNames)
output <- gjam(f$formula, xdata = f$xdata, ydata = f$ydata, 
                        modelList = ml)
#repeat with ng = 2000, burnin = 500, then:
pl  <- list(trueValues = f$trueValues, width = 3, height = 3)           
gjamPlot(output, pl)
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., D. Nemergut, B. Seyednasrollah, P. Turner, and S. Zhang. 2016.  Generalized joint attribute modeling for biodiversity analysis: Median-zero, multivariate, multifarious data.  Ecological Monographs, in press.
}
\seealso{
  \code{\link{gjamSimData}} simulates data
  \code{\link{gjam}} analyzes data
  
  A more detailed vignette is can be obtained with:

\code{browseVignettes('gjam')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

