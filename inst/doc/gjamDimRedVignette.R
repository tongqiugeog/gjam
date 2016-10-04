## ---- eval=FALSE---------------------------------------------------------
#  library(gjam)
#  S   <- 200
#  f   <- gjamSimData(n = 100, S = S, typeNames='CA')
#  rl  <- list(r = 5, N = 20)
#  ml  <- list(ng = 2000, burnin = 500, typeNames = f$typeNames,
#              reductList = rl, PREDICTX = F)
#  out <- gjamGibbs(f$formula, f$xdata, f$ydata, modelList = ml)
#  pl  <- list(trueValues = f$trueValues, SMALLPLOTS = F,
#              GRIDPLOTS=T, specLabs = F)
#  gjamPlot(output = out, plotPars = pl)

## ----fungal data summary, bty='n', fig.width=6.5, eval=F-----------------
#  library(gjam)
#  library(repmis)
#  source_data("https://github.com/jimclarkatduke/gjam/blob/master/fungEnd.RData?raw=True")
#  
#  xdata  <- fungEnd$xdata
#  otu    <- gjamReZero(fungEnd$yDeZero)
#  status <- fungEnd$status
#  
#  par(mfrow=c(1,3), bty='n', mar=c(1,1,1,1), oma = c(0,0,0,0),
#      mar = c(3,2,2,1), tcl = -0.5, mgp = c(3,1,0), family='')
#  hist(status, main='Host condition (morbid = 0)', ylab = 'Host obs')
#  hist(otu, nclass=100, ylab = 'Reads', main='each observation')
#  nobs <- gjamTrimY(otu, minObs = 1, OTHER = F)$nobs
#  hist(nobs, nclass=100, ylab = 'Total reads per OTU', main='Full sample')

## ----get y, eval = F-----------------------------------------------------
#  tmp <- gjamTrimY(otu, minObs = 100)
#  y   <- tmp$y
#  dim(fungEnd$y)               # all OTUs
#  dim(y)                       # trimmed data
#  tail(colnames(y))            # 'other' class added

## ----trim y, eval=FALSE--------------------------------------------------
#  ydata <- cbind(status, y) # host status is also a response
#  S     <- ncol(ydata)
#  typeNames    <- rep('CC',S)   # composition count data
#  typeNames[1] <- 'PA'          # binary host status

## ----factors, eval = F---------------------------------------------------
#  xdata$host <- relevel(xdata$host,'acerRubr')

## ----model fit, eval=F, eval=FALSE---------------------------------------
#  rl <- list(r = 5, N = 20)
#  ml <- list(ng = 2000, burnin = 500, typeNames = typeNames, reductList = rl)
#  output <- gjamGibbs(~ host*poly, xdata, ydata, modelList = ml)

## ----plots, eval=F-------------------------------------------------------
#  S <- ncol(ydata)
#  specColor     <- rep('black',S)
#  specColor[1]  <- 'red' # highlight host status
#  plotPars      <- list(corLines=F, specColor = specColor, GRIDPLOTS=T,
#                        specLabs = F, sdScaleY = T, SMALLPLOTS = F)
#  fit <- gjamPlot(output, plotPars)
#  fit$eComs[1:5,]

## ---- bstatus, eval=F----------------------------------------------------
#  beta <- fit$summaryCoeffs$betaCoeff
#  ws   <- grep('status_',rownames(beta))  # find coefficients for status
#  beta[ws,]

## ----bsig, eval=F--------------------------------------------------------
#  fit$summaryCoeffs$betaSig['status',]

## ----int, eval=F---------------------------------------------------------
#  xvector <- output$x[1,]*0
#  xnames  <- colnames(output$x)
#  names(xvector)  <- xnames
#  
#  xvector['hostfraxAmer'] <- 1
#  xvector['polypoly'] <- 1
#  fit1 <- gjamIIE(output, xvector, omitY = 'other')
#  
#  par(mfrow=c(1,3), bty='n', mar=c(1,1,1,1), oma = c(0,0,0,0),
#      mar = c(3,2,2,1), tcl = -0.5, mgp = c(3,1,0))
#  gjamIIEplot(fit1, response = 'status', effectMu = 'direct',
#              effectSd = 'direct',
#              legLoc = 'bottomright', ylim=c(-.5,.5))
#  title('Direct effect by host')
#  
#  gjamIIEplot(fit1, response = 'status', effectMu = 'int', effectSd = 'int',
#              legLoc = 'topright', ylim=c(-.5,.5))
#  title('Interactions with polyculture')
#  
#  gjamIIEplot(fit1, response = 'status', effectMu = 'ind', effectSd = 'ind',
#              legLoc = 'topright', ylim=c(-.5,.5))
#  title('Indirect effect of microbiome')

## ----predict, eval=F-----------------------------------------------------
#  y0 <- ydata[,1,drop=F]*0       #unhealthy host
#  
#  newdata   <- list(ydataCond = y0, nsim=50)
#  morbid    <- gjamPredict(output, newdata = newdata)
#  
#  newdata   <- list(ydataCond = y0 + 1, nsim = 50 )
#  healthy   <- gjamPredict(output, newdata = newdata)
#  
#  # compare predictions
#  par(mfrow=c(1,2), bty='n')
#  plot(healthy$sdList$yMu[,-1],morbid$sdList$yMu[,-1], cex=.4,
#       xlab='healty',ylab='morbid')
#  abline(0, 1, lty=2,col='grey')
#  plot(output$y[,2:20],healthy$sdList$yMu[,2:20], cex=.4,col='orange',
#       xlab='Observed',ylab='Predicted', pch=16)
#  points(output$y[,2:20],morbid$sdList$yMu[,2:20], cex=.4,col='blue', pch=16)
#  abline(0, 1, lty=2,col='grey')

