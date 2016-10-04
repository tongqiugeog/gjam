## ----traitBox1, fig.width = 7, fig.height = 3.5, echo = FALSE------------

.getBox <- function(x0,y0,tall,wide,col='white'){
  x1 <- x0 + wide
  y1 <- y0 + tall
  rect(x0, y0, x1, y1, col = col)
  mx <- mean(c(x0,x1))
  my <- mean(c(y0,y1))
  invisible(list( vec = c(x0,x1,y0,y1), mu = c(mx,my)) )
}

par(bty='n', mar=c(1,1,1,1), oma = c(0,0,0,0), 
    mar = c(3,2,2,1), tcl = -0.5, mgp = c(3,1,0), family='serif')

n <- 100
S <- 70
M <- 16
P <- 6
Q <- 14

xbox <- c(M,Q,M)
ybox <- c(n,n,Q)
xb <- c('M','Q','M')
yb <- c('n','n','Q')
xgap <- c(8,5,5)

ymax <- n + 6
xmax <- sum(xbox) + sum(xgap) + 5

plot(0,0,xlim=c(0,xmax),ylim=c(0,ymax),xaxt='n',yaxt='n',xlab='',ylab='',
     cex=.01)
xs <- 2

ti <- c('=','x','x')
ci <- c('U','X','beta')

for(j in 1:length(xbox)){
  
  ylo <- ymax - ybox[j]
  
  tmp <- .getBox(xs,ylo,ybox[j],xbox[j])
  xs  <- xgap[j] + tmp$vec[2]
  
  text(tmp$mu[1],ylo - 6,paste(yb[j],' x ',xb[j]))

  if(j < length(xbox))text(tmp$vec[2] + xgap[j]/2,ymax - ybox[j+1]/2, ti[j])
  if(j == 1)text(tmp$mu[1],tmp$mu[2],
                 expression(paste(italic(E),'[',bold(U),']')))
  if(j == 2)text(tmp$mu[1],tmp$mu[2],expression(bold(X)))
  if(j == 3)text(tmp$mu[1],tmp$mu[2],expression(hat(Alpha)))
}

## ----input, eval = T-----------------------------------------------------
library(gjam)
library(repmis)
source_data("https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True")

xdata <- forestTraits$xdata                    # n X Q
types <- forestTraits$traitTypes               # 12 trait types 
sbyt  <- forestTraits$specByTrait              # S X 12
pbys  <- gjamReZero(forestTraits$treesDeZero)  # n X S
head(sbyt)

## ----input2, eval = T----------------------------------------------------
table(sbyt$leaf)      # four levels

table(sbyt$xylem)     # diffuse/tracheid vs ring-porous

table(sbyt$repro)     # two levels

## ----input3, eval = T----------------------------------------------------
tmp         <- gjamSpec2Trait(pbys, sbyt, types)
tTypes      <- tmp$traitTypes                  # M = 15 values
u           <- tmp$plotByCWM                   # n X M
censor      <- tmp$censor                      # (0, 1) censoring, two-level CAT's
specByTrait <- tmp$specByTrait                 # S X M
M           <- ncol(u)
n           <- nrow(u)
types                                          # 12 individual trait types
cbind(colnames(u),tTypes)                      # M trait names and types

## ----setup2, eval = F----------------------------------------------------
#  censorList    <- gjamCensorY(values = c(0,1), intervals = cbind( c(-Inf,0),c(1,Inf) ),
#                               y = u, whichcol = c(13:14))$censor

## ----xdata, eval = T, echo = FALSE---------------------------------------
head(xdata)
table(xdata$soil)

## ----soil, eval = F------------------------------------------------------
#  xdata$soil <- relevel(xdata$soil,'reference')

## ----fit, eval = F-------------------------------------------------------
#  ml  <- list(ng = 3000, burnin = 500, typeNames = tTypes, holdoutN = 20,
#              censor=censor, notStandard = c('u1','u2','u3'))
#  out <- gjamGibbs(~ temp + stdage + moisture*deficit + deficit*soil,
#                   xdata = xdata, ydata = u, modelList = ml)
#  tnames    <- colnames(u)
#  specColor <- rep('black', M)                           # highlight types
#  wo <- which(tnames %in% c("leafN","leafP","SLA") )     # foliar traits
#  wf <- grep("leaf",tnames)                              # leaf habit
#  wc <- which(tnames %in% c("woodSG","diffuse","ring") ) # wood anatomy
#  
#  specColor[wc] <- 'brown'
#  specColor[wf] <- 'darkblue'
#  specColor[wo] <- 'darkgreen'
#  
#  pl  <- list(GRIDPLOTS = TRUE, plotAllY = T, specColor = specColor,
#              SMALLPLOTS = F, sigOnly=F, ncluster = 3)
#  fit <- gjamPlot(output = out, plotPars = pl)

## ----fit pars, eval = F--------------------------------------------------
#  out$modelSummary$betaMu      # Q by M coefficient matrix alpha
#  out$modelSummary$betaSe      # Q by M coefficient std errors
#  out$modelSummary$sigMu       # M by M covariance matrix omega
#  out$modelSummary$sigSe       # M by M covariance std errors

## ----fitTable, eval = F--------------------------------------------------
#  fit$betaEstimates[1:5,]      # Q by M coefficient matrix alpha

## ----IIEx, eval = F------------------------------------------------------
#  xdrydry <- xwetdry  <- out$x[1,]
#  xdrydry['moisture'] <- xdrydry['deficit'] <- -1
#  xwetdry['moisture'] <- 1
#  xwetdry['deficit']  <- -1

## ----IIE1, eval = F------------------------------------------------------
#  par(mfrow=c(2,2), bty='n', mar=c(1,3,1,1), oma = c(0,0,0,0),
#      mar = c(3,2,2,1), tcl = -0.5, mgp = c(3,1,0), family='')
#  
#  fit1 <- gjamIIE(output = out, xvector = xdrydry)
#  fit2 <- gjamIIE(output = out, xvector = xwetdry)
#  
#  gjamIIEplot(fit1, response = 'leafbroaddeciduous',
#              effectMu = c('main','int'),
#              effectSd = c('main','int'), legLoc = 'bottomleft',
#              ylim=c(-.31,.3)
#  title('deciduous')
#  gjamIIEplot(fit1, response = 'leafneedleevergreen',
#              effectMu = c('main','int'),
#              effectSd = c('main','int'), legLoc = 'bottomleft',
#              ylim=c(-.3,.3))
#  title('evergreen')
#  
#  gjamIIEplot(fit2, response = 'leafbroaddeciduous',
#              effectMu = c('main','int'),
#              effectSd = c('main','int'), legLoc = 'bottomleft',
#              ylim=c(-.3,.3))
#  gjamIIEplot(fit2, response = 'leafneedleevergreen',
#              effectMu = c('main','int'),
#              effectSd = c('main','int'), legLoc = 'bottomleft',
#              ylim=c(-.3,.3))

## ----IIE4, eval = F------------------------------------------------------
#  xvector <- out$x[1,]
#  par(mfrow=c(2,1), bty='n', mar=c(1,1,1,1), oma = c(0,0,0,0),
#      mar = c(3,2,2,1), tcl = -0.5, mgp = c(3,1,0), family='')
#  
#  omitY <- colnames(u)[colnames(u) != 'leafbroaddeciduous'] # omit all but deciduous
#  
#  fit <- gjamIIE(out, xvector)
#  gjamIIEplot(fit, response = 'leafP', effectMu = c('main','ind'),
#              effectSd = c('main','ind'), legLoc = 'topright',
#              ylim=c(-.6,.6))
#  title('foliar P')
#  gjamIIEplot(fit, response = 'leafN', effectMu = c('main','ind'),
#              effectSd = c('main','ind'), legLoc = 'bottomright',
#              ylim=c(-.6,.6))
#  title('foliar N')

## ----traitBox2, fig.width = 6.7, fig.height = 4, echo = FALSE------------
par(bty='n', bty='n', mar=c(1,1,1,1), oma = c(0,0,0,0), 
    mar = c(3,2,2,1), tcl = -0.5, mgp = c(3,1,0), family='serif')

n <- 100
S <- 70
M <- 16
P <- 6
Q <- 14

xbox <- c(S,M,Q,S,M)
ybox <- c(n,S,n,Q,S)
xb   <- c('S','M','Q','S','M')
yb   <- c('n','S','n','Q','S')
xgap <- c(15,28,15,15,15,15)

ymax <- n + 5
xmax <- sum(xbox) + sum(xgap) + 5

plot(0,0,xlim=c(0,xmax),ylim=c(0,ymax),xaxt='n',yaxt='n',xlab='',ylab='',
     cex=.01)
xs <- xgap[1]

ti <- c('x','=','x','x')
ci <- c('W','T','X','beta','T')
col <- rep('white',length(xbox))
col[c(2,5)] <- 'wheat'

for(j in 1:length(xbox)){
  
  ylo <- ymax - ybox[j]
  tmp <- .getBox(xs,ylo,ybox[j],xbox[j], col[j])
  xs  <- xgap[j] + tmp$vec[2]
  
  text(tmp$mu[1], ylo - 6, paste(yb[j],' x ',xb[j]) )

  if(j < length(xbox))text(tmp$vec[2] + xgap[j]/2,ymax - ybox[j+1]/2, ti[j])
  if(j == 1)text(tmp$mu[1],tmp$mu[2],
                 expression(paste(italic(E),'[',bold(W),']')))
  if(j == 2)text(tmp$mu[1],tmp$mu[2],expression(bold(T)))
  if(j == 3)text(tmp$mu[1],tmp$mu[2],expression(bold(X)))
  if(j == 4)text(tmp$mu[1],tmp$mu[2],expression(hat(beta)))
  if(j == 5)text(tmp$mu[1],tmp$mu[2],expression(bold(T)))
}

## ----PTM, eval = F-------------------------------------------------------
#  tl  <- list(plotByTrait = u, traitTypes = tTypes, specByTrait = specByTrait)
#  rl  <- list(r = 8, N = 20)
#  ml  <- list(ng = 1000, burnin = 200, typeNames = 'CC', holdoutN = 20,
#                    traitList = tl, reductList = rl)
#  out <- gjamGibbs(~ temp + stdage + deficit*soil, xdata = xdata,
#                       ydata = pbys, modelList = ml)
#  S <- nrow(specByTrait)
#  specColor <- rep('black',S)
#  
#  wr <- which(specByTrait[,'ring'] == 1)                  # ring porous
#  wb <- which(specByTrait[,'leafneedleevergreen'] == 1)   # evergreen
#  ws <- which(specByTrait[,'shade'] >= 4)                 # shade tolerant
#  specColor[wr] <- 'brown'
#  specColor[ws] <- 'black'
#  specColor[wb] <- 'darkgreen'
#  
#  par(family = '')
#  pl  <- list(width=4, height=4, corLines=F, SMALLPLOTS=F,GRIDPLOTS=T,
#                    specColor = specColor, ncluster = 8)
#  fit <- gjamPlot(output = out, pl)

## ----trait pars, eval = F------------------------------------------------
#  out$modelSummary$betaTraitMu   # Q by M coefficient matrix alpha
#  out$modelSummary$betaTraitSe   # Q by M coefficient std errors
#  out$modelSummary$sigmaTraitMu  # M by M covariance matrix omega
#  out$modelSummary$sigmaTraitSe  # M by M covariance std errors

## ----trait pred, eval = F------------------------------------------------
#  out$modelSummary$tMu[1:5,]     # n by M predictive means
#  out$modelSummary$tSd[1:5,]     # n by M predictive std errors

## ----ecoms, eval = F-----------------------------------------------------
#  fit$eComs[,1:4]

