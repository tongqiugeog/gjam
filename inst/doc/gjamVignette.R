## ----fig1, fig.width = 6.7, fig.height = 4, echo = FALSE-----------------

  sig    <- .9
  mu     <- 3.1
  offset <- -2
  
  par(mfrow = c(1, 2), bty = 'n', mar = c(4, 5, 3, .1), cex=1.3, family='serif')
  part <- c(0, 2.2, 3.3, 4.5, 6.6)
  w    <- seq(-1, 7, length = 1000)
  dw   <- dnorm(w, mu, sig)
  dp   <- dw[ findInterval(part, w) ]
  pw   <- pnorm(part, mu, sig)
  pw[-1] <- diff(pw)
  
  plot(w, 2*dw - .5, type = 'l', ylim = c(-.5, 4), yaxt = 'n', 
       ylab = expression(paste(italic(y), '|', italic(w), ', ', bold(p), 
                               sep = '')), 
       xlab = expression(paste(italic(w), '|', bold(x), ', ', bold(beta), 
                              ', ', bold(Sigma), sep = '')), 
       xlim = c(offset, 7), lwd = 2)
  axis(2, at = c(0:5))
  
  db <- .15
  int <- 4
  
  polygon( c(w, rev(w)), 2*c(dw, w*0) - .5, col = 'grey', lwd = 2)
  lines(c(-1, part[1]), c(0, 0), lwd = 2)
  
  for(j in 1:(length(part))){
    
    lines( part[j:(j+1)], c(j, j), lwd = 3)
    ww <- which(w >= part[j] & w <= part[j+1])
    
    if(j == 3){
      w1 <- ww[1]
      w2 <- max(ww)
      arrows( mean(w[ww]), 2*max(dw[ww]) - .4, mean(w[ww]), 
              j - .4, angle = 20, lwd = 5, col = 'grey', length = .2)
      arrows( w[w1] - .5 , j , -.7, j , angle = 20, 
              lwd = 5, col = 'grey', length = .2)
      text( c(w[w1], w[w2]), c(3.3, 3.3), 
            expression(italic(p)[4], italic(p)[5]), cex=.9)
      text( w[w2] + .3, .6, expression( italic(w)[italic(is)] ))
      text( 0, 3.5, expression( italic(y)[italic(is)] ))
    }
    
    coll <- 'white'
    if(j == int)coll <- 'grey'
    rect( offset, j - 1 - db, 2*pw[j] + offset, j - 1 + db, 
          col = coll, border = 'black', lwd = 2)
  }
  
  ww <- which(w >= part[int - 1] & w <= part[int])
  abline(h = -.5, lwd = 2)
  
  title('a) Data generation', adj = 0, font.main = 1, font.lab = 1, cex=.8)
  
  plot(w, 2*dw - .5, type = 'l', ylim = c(-.5, 4), yaxt = 'n', 
       ylab = expression(italic(y)), 
       xlab = expression(paste(italic(w), '|', italic(y), ', ', bold(p), sep = '')), 
       xlim = c(offset, 7), lwd = 2, col = 'grey')
  axis(2, at = c(0:5))
  
  abline(h = -.5, lwd = 2, col = 'grey')
  
  wseq <- c(-10,part)
  for(j in 1:(length(part))){
    
    coll <- 'white'
    border <- 'grey'
    
    if(j == int){
      coll <- 'grey'
      border <- 'black'
      rect( offset, j - 1 - db, 2*pw[j] + offset, j - 1 + db, 
            col = 'black', border = 'black')
    }
    lines( part[j:(j+1)], c(j, j), lwd = 3)
    lines(part[c(j, j)], 2*c(0, dp[j])-.5, col = 'grey')
  }
  
  lines(c(-1, part[1]), c(0, 0), lwd = 2)
  ww <- which(w >= part[int - 1] & w <= part[int])
  polygon( w[c(ww, rev(ww))], 2*c(dw[ww], ww*0) - .5, col = 'grey', lwd = 2)
  
  arrows( mean(w[ww]),  int - 1.3, mean(w[ww]),  2*max(dw) - .5, 
          angle = 20, lwd = 5, col = 'grey', length = .2)
  arrows( -.5,  int - 1, min(w[ww]) - .4, int - 1, angle = 20, 
          lwd = 5, col = 'grey', length = .2)
  
  title('b) Inference', adj = 0, font.main = 1, font.lab = 1, cex=.8)

## ----simulate CA, eval = F-----------------------------------------------
#  library(gjam)
#  f <- gjamSimData(n = 500, S = 10, Q = 3, typeNames = 'CA')
#  summary(f)

## ----show formula, eval = F----------------------------------------------
#  f$formula

## ----plot simulated y, fig.show = "hold", fig.width = 6.5, eval = F------
#  par(bty = 'n', mfrow = c(1,2), family='')
#  h <- hist(c(-1,f$y),nclass = 50,plot = F)
#  plot(h$counts,h$mids,type = 's')
#  plot(f$w,f$y,cex = .2)

## ----fit CA data, eval = F-----------------------------------------------
#  # a few iterations
#  ml  <- list(ng = 1000, burnin = 100, typeNames = f$typeNames)
#  out <- gjamGibbs(f$formula, f$xdata, f$ydata, modelList = ml)
#  summary(out)

## ----summary of chains, eval = F-----------------------------------------
#  summary(out$chains)

## ----summary of fitted model, eval = FALSE-------------------------------
#  summary(out$modelSummary)

## ----show classes, eval = F----------------------------------------------
#  out$modelSummary$classBySpec

## ----betaMu, eval=F------------------------------------------------------
#  out$parameterTables$betaMu

## ----betaSd, eval=F------------------------------------------------------
#  out$parameterTables$betaSe

## ----plot CA data, family='', eval = FALSE-------------------------------
#  f   <- gjamSimData(n = 500, S = 10, typeNames = 'CA')
#  ml  <- list(ng = 1000, burnin = 200, typeNames = f$typeNames)
#  out <- gjamGibbs(f$formula, f$xdata, f$ydata, modelList = ml)
#  pl  <- list(trueValues = f$trueValues, GRIDPLOTS = T, SMALLPLOTS = F)
#  gjamPlot(output = out, plotPars = pl)

## ----example output, fig.show = "hold", fig.width = 6.5, eval = F--------
#  par(bty = 'n', mfrow = c(1,3), family='')
#  plot(f$trueValues$beta, out$parameterTables$betaMu, cex = .2)
#  plot(f$trueValues$corSpec, out$parameterTables$corMu, cex = .2)
#  plot(f$y,out$modelSummary$ypredMu, cex = .2)

## ----design1, eval = F---------------------------------------------------
#  library(repmis)
#  d <- "https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True"
#  source_data(d)
#  xdata <- forestTraits$xdata[,c(1,2,8)]

## ----setupX, echo=F, eval=T----------------------------------------------
temp <- c(1.22,  0.18, -0.94,  0.64,  0.82)
deficit <- c(0.04,  0.21,  0.20,  0.82, -0.18)
soil <- c('reference', 'reference', 'SpodHist',  'reference', 'reference')
xx <- data.frame( temp, deficit, soil )
attr(xx$soil,'levels') <- c("reference","SpodHist","EntVert","Mol","UltKan")

## ----design1.2, eval = F-------------------------------------------------
#  xdata[1:5,]

## ----design2, eval = T---------------------------------------------------
formula <- as.formula( ~ temp + deficit + soil )

## ----design3, echo=F-----------------------------------------------------
  tmp <- model.frame(formula,data=xx)
  x   <- model.matrix(formula, data=tmp)
  x[1:5,]

## ----design4, eval = T---------------------------------------------------
formula <- as.formula( ~ temp*soil )

## ----design5, echo = F---------------------------------------------------
tmp <- model.frame(formula,data=xx,na.action=NULL)
x   <- model.matrix(formula, data=tmp)
x[1:5,]

## ----design6, eval = T---------------------------------------------------
formula <- as.formula( ~ temp + I(temp^2) + deficit )

## ----design7, echo = F---------------------------------------------------
tmp <- model.frame(formula,data=xx,na.action=NULL)
x   <- model.matrix(formula, data=tmp)
x[1:5,]

## ----design8, eval = T---------------------------------------------------
formula <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )

## ----design9, echo = F---------------------------------------------------
tmp <- model.frame(formula,data=xx,na.action=NULL)
x   <- model.matrix(formula, data=tmp)
  x[1:5,]

## ----get trees, eval = F-------------------------------------------------
#  library(gjam)
#  ydata  <- gjamReZero(forestTraits$treesDeZero)  # extract y
#  dim(ydata)
#  ydata[1:5,1:6]

## ----fit trees, eval = F-------------------------------------------------
#  rl   <- list(r = 8, N = 20)
#  ml   <- list(ng = 1000, burnin = 100, typeNames = 'DA', reductList = rl)
#  form <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )
#  out  <- gjamGibbs(form, xdata = xdata, ydata = ydata, modelList = ml)
#  pl   <- list(SMALLPLOTS = F, GRIDPLOTS=T, corLines=F, specLabs = F)
#  gjamPlot(output = out, plotPars = pl)

## ----plot save, eval = F-------------------------------------------------
#  plotPars <- list(SMALLPLOTS = F, GRIDPLOTS=T, SAVEPLOTS = T, outfolder = 'plots')

## ----effort simulation, eval = F-----------------------------------------
#  S   <- 5
#  n   <- 50
#  ef  <- list( columns = 1:S, values = round(runif(n,.5,5),1) )
#  f   <- gjamSimData(n, S, typeNames = 'DA', effort = ef)
#  ef

## ----w vs y, bty = 'n', fig.width = 6.5, eval = F------------------------
#  plot(f$w,f$y, cex = .2)

## ----including effort, bty = 'n', fig.width = 6.5, eval = F--------------
#  plot(f$w*ef$values, f$y, cex = .2)

## ----fitting, eval = F---------------------------------------------------
#  S   <- 10
#  n   <- 1500
#  ef  <- list( columns = 1:S, values = round(runif(n,.5,5),1) )
#  f   <- gjamSimData(n, S, typeNames = 'DA', effort = ef)
#  ml  <- list(ng = 1000, burnin = 250, typeNames = f$typeNames, effort = ef)
#  out <- gjamGibbs(f$formula, f$xdata, f$ydata, modelList = ml)
#  pl  <- list(trueValues = f$trueValues,SMALLPLOTS=F)
#  gjamPlot(output = out, plotPars = pl)

## ----betaPrior, eval = F-------------------------------------------------
#  source_data("https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True")
#  
#  xdata <- forestTraits$xdata
#  y     <- gjamReZero(forestTraits$treesDeZero)
#  ydata <- gjamTrimY(y,300)$y        # a sample of species
#  types <- 'DA'
#  
#  xnames <- c('temp','deficit')      # variables for truncated priors
#  Q      <- length(xnames)
#  S      <- ncol(ydata)
#  
#  loBeta <- matrix(-Inf,Q,S)         # initialize priors
#  hiBeta <- matrix(Inf,Q,S)
#  rownames(loBeta) <- rownames(hiBeta) <- xnames
#  
#  loBeta['temp',]    <- 0            # minimum zero
#  hiBeta['deficit',] <- 0            # maximum zero
#  
#  bp <- list(lo = loBeta, hi = hiBeta)
#  rl <- list(N = 10, r = 5)          # dimension reduction
#  modelList <- list(ng = 5000, burnin = 500, typeNames = types,
#                    betaPrior = bp, reductList = rl)

## ----compData, eval = FALSE----------------------------------------------
#  f     <- gjamSimData(S = 8, typeNames = 'CC')
#  ml    <- list(ng = 2000, burnin = 500, typeNames = f$typeNames)
#  out   <- gjamGibbs(f$formula, f$xdata, f$ydata, modelList = ml)
#  pl    <- list(trueValues = f$trueValues, SMALLPLOTS = F)
#  gjamPlot(output = out, plotPars = pl)

## ----compFC, eval = FALSE------------------------------------------------
#  f     <- gjamSimData(S = 20, typeNames = 'FC')
#  ml    <- list(ng = 2000, burnin = 500, typeNames = f$typeNames)
#  out   <- gjamGibbs(f$formula, f$xdata, f$ydata, modelList = ml)
#  pl    <- list(trueValues = f$trueValues, SMALLPLOTS = F)
#  gjamPlot(output = out, plotPars = pl)

## ----ordinal, eval = FALSE-----------------------------------------------
#  f   <- gjamSimData(typeNames = 'OC')
#  ml  <- list(ng = 2000, burnin = 500, typeNames = f$typeNames)
#  out <- gjamGibbs(f$formula, f$xdata, f$ydata, modelList = ml)
#  print(out)

## ----ordinal partition, eval = FALSE-------------------------------------
#  keep <- strsplit(colnames(out$parameterTables$cutMu),'C-') #only saved columns
#  keep <- matrix(as.numeric(unlist(keep)), ncol = 2, byrow = T)[,2]
#  plot(f$trueValues$cuts[,keep],out$parameterTables$cutMu)

## ----ordPlots, eval = FALSE----------------------------------------------
#  pl  <- list(trueValues = f$trueValues, SMALLPLOTS = F)
#  gjamPlot(output = out, plotPars = pl)

## ----cat, eval = T, echo = F---------------------------------------------
leaf  <- c('bd', 'nd', 'be', 'other')
xylem <- c('dp', 'rp', 'other') 
np    <- 7
xx <- data.frame( leaf = factor(sample(leaf,np,replace=T)),
            xylem = factor(sample(xylem,np,replace=T) ))
xx

## ----catY, eval = T, echo = F--------------------------------------------
wl <- match(xx$leaf,leaf)
wx <- match(xx$xylem,xylem)
ml <- matrix(0,np,4)
ml[cbind(1:np,wl)] <- 1
colnames(ml) <- paste('leaf',leaf,sep='_')
mx <- matrix(0,np,3)
mx[cbind(1:np,wx)] <- 1
colnames(mx) <- paste('xylem',xylem,sep='_')
cbind(ml,mx)

## ----cat2, eval = FALSE--------------------------------------------------
#  types <- c('CAT','CAT')
#  f     <- gjamSimData(n=2000, S = length(types), typeNames = types)
#  ml    <- list(ng = 1500, burnin = 500, typeNames = f$typeNames, PREDICTX = F)
#  out   <- gjamGibbs( f$formula, xdata = f$xdata, ydata = f$ydata, modelList = ml )
#  pl    <- list(trueValues = f$trueValues, SMALLPLOTS=F, plotAllY = T)
#  gjamPlot(out, plotPars = pl)

## ----many types, eval = FALSE--------------------------------------------
#  types <- c('OC','OC','OC','OC','CC','CC','CC','CC','CC','CA','CA','PA','PA')
#  f     <- gjamSimData(S = length(types), Q = 3, typeNames = types)
#  ml    <- list(ng = 2000, burnin = 500, typeNames = f$typeNames)
#  out   <- gjamGibbs(f$formula, f$xdata, f$ydata, modelList = ml)
#  tmp   <- data.frame(f$typeNames, out$modelSummary$classBySpec[,1:10])
#  print(tmp)

## ----mixed analysis, eval = FALSE----------------------------------------
#  pl  <- list(trueValues = f$trueValues, SMALLPLOTS = F)
#  gjamPlot(output = out, plotPars = pl)

## ----simulate missing data, eval = FALSE---------------------------------
#  f <- gjamSimData(typeNames = 'OC', nmiss = 20)
#  which(is.na(f$xdata), arr.ind = T)

## ----holdouts, eval = FALSE----------------------------------------------
#  f   <- gjamSimData(typeNames = 'CA', nmiss = 20)
#  ml  <- list(ng = 2000, burnin = 500, typeNames = f$typeNames, holdoutN = 50)
#  out <- gjamGibbs(f$formula, f$xdata, f$ydata, modelList = ml)
#  
#  par(mfrow=c(1,2))
#  xMu  <- out$modelSummary$xpredMu
#  xSd  <- out$modelSummary$xpredSd
#  yMu  <- out$modelSummary$ypredMu
#  hold <- out$holdoutIndex
#  
#  plot(out$x[hold,-1],xMu[hold,-1], cex=.2)
#  title('holdouts in x'); abline(0,1)
#  plot(out$y[hold,], yMu[hold,], cex=.2)
#  title('holdouts in y'); abline(0,1)

## ----effortPredict, eval = FALSE-----------------------------------------
#  sc  <- 3                               #no. CA responses
#  sd  <- 10                              #no. DA responses
#  tn  <- c( rep('CA',sc),rep('DA',sd) )  #combine CA and DA obs
#  S   <- length(tn)
#  n   <- 500
#  emat   <- matrix( runif(n,.5,5), n, sd)              #simulated DA effort
#  effort <- list(columns = c((sc+1):S), values = emat )
#  f      <- gjamSimData(n = n, typeNames = tn, effort = effort)
#  ml     <- list(ng = 2000, burnin = 500, typeNames = f$typeNames, effort = f$effort)
#  out    <- gjamGibbs(f$formula, f$xdata, f$ydata, modelList = ml)
#  
#  par(mfrow=c(1,2),bty='n')
#  gjamPredict(out, y2plot = colnames(f$ydata)[tn == 'DA']) #predict DA data

## ----effortPredictNew, eval = FALSE--------------------------------------
#  newdata   <- list(xdata = f$xdata, effort=effort, nsim = 50 )      # effort unchanged
#  p1 <- gjamPredict(output = out, newdata = newdata)
#  
#  plot(f$y[,tn == 'DA'], p1$sdList$yMu[,tn == 'DA'],ylab = 'Predicted',cex=.1)
#  abline(0,1)
#  
#  newdata$effort$values <- effort$values*0 + 1       # predict for effort = 1
#  p2 <- gjamPredict(output = out, newdata = newdata)
#  
#  points(f$y[,tn == 'DA'], p2$sdList$yMu[,tn == 'DA'],col='orange',cex=.1)
#  abline(0,1)

## ----effortPredictCond, eval = FALSE-------------------------------------
#  newdata <- list(ydataCond = f$y[,1:2], nsim=200)   # cond on obs CA data
#  p1      <- gjamPredict(output = out, newdata = newdata)$sdList$yMu[,tn == 'DA']
#  
#  yc     <- f$y[,1:2]                                  # cond on new CA values
#  yc[,1] <- mean(yc[,1])
#  yc[,2] <- 0
#  newdata   <- list(ydataCond = yc, nsim=200)
#  p2 <- gjamPredict(output = out, newdata = newdata)$sdList$yMu[,tn == 'DA']
#  plot(f$y[,tn == 'DA'], p1, xlab='Obs', ylab = 'Pred', cex=.1, ylim=range(c(p1,p2)))
#  points(f$y[,tn == 'DA'], p2,col='orange',cex=.1)
#  abline(0,1)

## ----cont1, echo=F-------------------------------------------------------
D <- rbind( c(1, -1, -1), c(0,1,0), c(0, 0, 1))
colnames(D) <- c('intercept','b','c')
rownames(D) <- c('a','b','c')
C <- D
C[,1] <- 1
C

## ----cont2, echo=F-------------------------------------------------------
t(D)

