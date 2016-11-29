
.combineFacLevels <- function(xfactor,fname=NULL, aname = 'reference', 
                              vminF=1){
  tmp <- as.character(xfactor)
  tmp[tmp %in% fname] <- aname
  tab <- table(tmp)
  wm  <- names(tab)[tab < vminF]
  tmp[tmp %in% wm] <- aname
  as.factor(tmp)
}

.factor2Numeric <- function(xfactor) as.numeric(as.character(xfactor))

.getColor <- function(col,trans){
  
  # trans - transparency fraction [0, 1]
  
  tmp <- col2rgb(col)
  rgb(tmp[1,], tmp[2,], tmp[3,], maxColorValue = 255, 
      alpha = 255*trans, names = paste(col,trans,sep='_'))
}
  
.figure1 <- function(){
  sig <- .9
  mu  <- 3.1
  offset <- -2
  
  par(mfrow=c(1,2),bty='n',mar=c(6,5,3,.1))
  part <- c(0,2,3.3,4.9,6.6)
  w    <- seq(-1,7,length=500)
  dw   <- dnorm(w,mu,sig)
  dp   <- dw[ findInterval(part,w) ]
  pw   <- pnorm(part,mu,sig)
  pw[-1] <- diff(pw)
  
  plot(w,2*dw - .5,type='l',ylim=c(-.5,4),yaxt='n', 
       ylab= expression(paste(italic(y),'|(',italic(w),', ',bold(p),')',sep='')), 
       xlab= expression(paste(italic(w),'|(',bold(x)[i],', ',bold(beta),
                              ', ',bold(Sigma),')',sep='')), 
       xlim=c(offset,7), lwd=2)
  axis(2, at = c(0:5))
  
  db <- .15
  int <- 4
  
  polygon( c(w,rev(w)), 2*c(dw,w*0) - .5, col='grey', lwd=2)
  
  lines(c(-1,part[1]),c(0,0),lwd=2)
  
  for(j in 1:(length(part))){
    
    lines( part[j:(j+1)],c(j,j), lwd=3)
    
    ww <- which(w >= part[j] & w <= part[j+1])
    
    if(j == 3){
      w1 <- ww[1]
      w2 <- max(ww)
      arrows( mean(w[ww]), 2*max(dw[ww]) - .4, mean(w[ww]), 
              j - .4, angle=20,lwd=3, col = 'grey', length=.2)
      arrows( w[w1] - .5 , j , -.7, j , angle= 20, 
              lwd = 3, col='grey', length=.2)
      text( c(w[w1], w[w2]),c(3.3,3.3),
            expression(italic(p)[4], italic(p)[5]))
      text( w[w2] + .3,.6,expression( italic(w)[italic(is)] ))
      text( 0,3.5,expression( italic(y)[italic(is)] ))
    }
    
    coll <- 'white'
    if(j == int)coll <- 'grey'
    rect( offset, j - 1 - db, 2*pw[j] + offset, j - 1 + db, 
          col=coll, border='black', lwd=2)
  }
  
  ww <- which(w >= part[int - 1] & w <= part[int])
  abline(h = -.5, lwd = 2)
  
  title('a) Data generation',adj=0, font.main = 1, font.lab =1)
  
  plot(w,2*dw - .5,type='l',ylim=c(-.5,4), yaxt='n', 
       ylab= expression(italic(y)), 
       xlab= expression(paste(italic(w),'|(',italic(y),', ',bold(p),')',sep='')), 
       xlim=c(offset,7), lwd=2,col='grey')
  axis(2, at = c(0:5))
  
  lines(c(-1,part[1]),c(0,0),lwd=2)
  abline(h=-.5, lwd=2, col='grey')
  
  for(j in 1:(length(part))){
    
    lines( part[j:(j+1)],c(j,j), lwd=3)
    lines(part[c(j,j)],2*c(0,dp[j])-.5, col='grey')
    
    coll <- 'white'
    if(j -- int)coll <- 'grey'
    
    if(j == int){
      rect( offset, j - 1 - db, 2*pw[j] + offset, j - 1 + db,
            col='black', border='black')
    }
  }
  
  ww <- which(w >= part[int - 1] & w <= part[int])
  polygon( w[c(ww,rev(ww))], 2*c(dw[ww],ww*0) - .5, col='grey', lwd=2)
  
  arrows( mean(w[ww]),  int - 1.3,mean(w[ww]),  2*max(dw) - .5,
          angle=20,lwd=3, col = 'grey', length=.2)
  arrows( -.5,  int - 1, min(w[ww]) - .4, int - 1, angle= 20,
          lwd = 3, col='grey', length=.2)
  
  title('b) Inference',adj=0, font.main = 1, font.lab = 1)
}

.add2matrix <- function(values,xmat=NULL){
  
  #xmat   - n X ? matrix with one row, columns are integer values
  #values - length-n vector be added/slotted in to xvec
  
  if(is.null(xmat)){
    n    <- length(values)
    cc   <- sort(unique(values))
    xmat <- matrix(0,n,length(cc),dimnames = list(1:n,cc))
    xmat[ cbind( c(1:n),match(values,cc)) ] <- 1
    return(xmat)
  }
  
  n <- nrow(xmat)
  if(length(values) != n)stop('vector length must equal rows in xmat')
  
  all <- sort( unique( c(values,as.numeric(colnames(xmat))) ))
  nc       <- length(all)
  
  xnew <- matrix(0,n,nc,dimnames = list(1:n,all))
  xnew[,colnames(xmat)] <- xmat
  
  xnew[ cbind(c(1:n),match(values,all)) ] <- xnew[ cbind(c(1:n),match(values,all)) ] + 1
  xnew
}
.appendMatrix <- function(m1,m2,fill=NA,SORT=F,asNumbers=F){  
  
  # matches matrices by column names
  # asNumbers: if column heads are numbers and SORT, then sort numerically

   if(length(m1) == 0){
     if(is.matrix(m2)){
       m3 <- m2
     } else {
       m3 <- matrix(m2,nrow=1)
     }
     if( !is.null(names(m2)) )colnames(m3) <- names(m2)
     return(m3)
   }
   if(length(m2) == 0){
     if(!is.matrix(m1))m1 <- matrix(m1,nrow=1)
     return(m1)
   }
   if( is.vector(m1) | (length(m1) > 0 & !is.matrix(m1)) ){
     nn <- names(m1)
     if(is.null(nn))warning('cannot append matrix without names')
     m1 <- matrix(m1,1)
     colnames(m1) <- nn
   }  
   if( is.vector(m2) | (length(m2) > 0 & !is.matrix(m2)) ){
     nn <- names(m2)
     if(is.null(nn))warning('cannot append matrix without names')
     m2 <- matrix(m2,1)
     colnames(m2) <- nn
   }

   c1 <- colnames(m1)
   c2 <- colnames(m2)
   r1 <- rownames(m1)
   r2 <- rownames(m2)
   n1 <- nrow(m1)
   n2 <- nrow(m2)

   allc <-  unique( c(c1,c2) ) 
   if(SORT & !asNumbers)allc <- sort(allc)
   if(SORT & asNumbers){
     ac <- as.numeric(allc)
     allc <- as.character( sort(ac) )
   }

   nr <- n1 + n2
   nc <- length(allc)

   if(is.null(r1))r1 <- paste('r',c(1:n1),sep='-')
   if(is.null(r2))r2 <- paste('r',c((n1+1):nr),sep='-')
   new <- c(r1,r2)

   mat1 <- match(c1,allc)
   mat2 <- match(c2,allc)

   out <- matrix(fill,nr,nc)
   colnames(out) <- allc
   rownames(out) <- new

   out[1:n1,mat1] <- m1
   out[(n1+1):nr,mat2] <- m2
   out
}
.byIndex <- function(xx,INDICES,FUN,coerce=F,...){  
  
#INDICES is list, each same length as  x
  
#  fun <- match.fun(FUN)
  
  nl <- length(INDICES)
  
  tmp  <-  unlist(by( as.vector(xx),INDICES,FUN,...) ) 
  nd   <- dim(tmp)
  tmp  <- array(tmp,dim=nd, dimnames=dimnames(tmp))
  
  tmp[is.na(tmp)] <- 0
  
  if(!coerce)return(tmp)
  
  dname <- dimnames(tmp)
  mk    <- rep(0,length(nd))
  
  for(k in 1:length(nd))mk[k] <- max(as.numeric(dimnames(tmp)[[k]]))
  
  wk <- which(mk > nd)
  if(length(wk) > 0){
    tnew  <- array(0,dim=mk)
    if(length(dim(tnew)) == 1)tnew <- matrix(tnew,dim(tnew),1)
    for(k in wk){
      newk <- c(1:mk[k])
      mat  <- match(dimnames(tmp)[[k]],newk)
      if(k == 1){
        tnew[mat,] <- tmp
        rownames(tnew) <- 1:nrow(tnew)
      }
      if(k == 2){
        tnew[,mat] <- tmp
        colnames(tnew) <- c(1:ncol(tnew))
      }
      tmp <- tnew
    }
  }
  tmp
}

.chains2density <- function(chainMat,labs=NULL,reverseM=F,varName=NULL,
                            cut=0){
  
  #assumes column names are varName or 'something_varname'
  
  #chainMat - MCMC output [samples,chains]
  
  chNames <- colnames(chainMat)
  
  if(!is.null(varName)){
    
    wc <- grep(varName,colnames(chainMat),fixed=T)
    if(length(wc) == 0)stop('varName not found in colnames(chainMat)')
    
    ww <- grep('_',colnames(chainMat),fixed=T)
    if(length(ww) > 0){
      tmp <- matrix( unlist(strsplit(colnames(chainMat),'_')),ncol=2,byrow=T) 
      wc <- which(tmp[,2] == varName)
      if(length(wc) == 0)wc <- which(tmp[,1] == varName)
    }
    chainMat <- chainMat[,wc]
    if(!is.matrix(chainMat))chainMat <- matrix(chainMat,ncol=1)
    colnames(chainMat) <- chNames[wc]
  }
  
  nj <- ncol(chainMat)
  nd <- 512
  
  clab <- colnames(chainMat)
  if(is.null(labs) & !is.null(clab))labs <- clab
  if(is.null(labs) & is.null(clab)) labs <- paste('v',c(1:nj),sep='-')
  
  xt <- yt <- matrix(NA,nj,nd)
  rownames(xt) <- rownames(yt) <- labs
  
  xrange <- signif(range(chainMat),2)
  
  for(j in 1:nj){
    
 #   lj  <- labs[j]
    xj  <- chainMat[,j]
    tmp <- density(xj,n = nd, cut=cut, na.rm=T)
    xt[j,]  <- tmp$x
    yt[j,]  <- tmp$y
    
  }
  yymax <- max(yt,na.rm=T)
  
  if(reverseM){
    xt <- -t( apply(xt,1,rev) )
    yt <- t( apply(yt,1,rev) )
  }
  
  list(x = xt, y = yt, xrange = xrange, ymax = yymax, chainMat = chainMat)
}
.checkDesign <- function(x,intName='intercept',xflag=':', 
                         isFactor = character(0)){  # 

  # xflag - indicates that variable is an interaction
  # isFactor - character vector of factor names returned if not supplied
  
  p <- ncol(x)
  
  if(ncol(x) < 3){
 #   message('no design check, < 2 covariates in x')
    return( list(VIF = 0, correlation = 1, rank = 2, p = 2, isFactor=isFactor) )
  }
    
  if(is.null(colnames(x))){
    colnames(x) <- paste('x',c(1:p),sep='_')
  }
  xrange      <- apply(x,2,range,na.rm=T)
  wi          <- which(xrange[1,] == 1 & xrange[2,] == 1)
  if(length(wi) > 0)colnames(x)[wi] <- 'intercept'
  
  wx <- grep(xflag,colnames(x))
  wi <- which(colnames(x) == 'intercept')
  wi <- unique(c(wi,wx))

  xname <- colnames(x)
  
  wmiss <- which(is.na(x),arr.ind=T)
  
  if(length(wmiss) > 0){
    rowTab <- table( table(wmiss[,1]) )
    colTab <- table(wmiss[,2])
  }
    
  VIF <- rep(NA,p)
  names(VIF) <- xname
  
  GETF <- F
  if(length(isFactor) > 0)GETF <- T
  
  for(k in 1:p){

    if(xname[k] %in% wi)next
    
    notk <- xname[xname != xname[k] & !xname %in% xname[wi]]
    ykk  <- x[,xname[k]]
    xkk  <- x[,notk,drop=F]
    
    wna <- which(is.na(ykk) | is.na(rowSums(xkk)))
    if(length(wna) > 0){
      ykk <- ykk[-wna]
      xkk <- xkk[-wna,]
    }
    
    tkk <- summary(lm(ykk ~ xkk))$adj.r.squared
    VIF[k] <- 1/(1 - tkk)
    
    xu <- sort( unique(x[,k]) )
    tmp <- identical(c(0,1),xu)
    if(GETF)if(tmp)isFactor <- c(isFactor,xname[k])
  }

  VIF <- VIF[-wi] 

  corx <- cor(x[,-wi], use="complete.obs")
  if(length(wna) == 0){
    rankx <- qr(x)$rank
  } else {
    rankx <- qr(x[-wna,])$rank
  }
  corx[upper.tri(corx,diag=T)] <- NA
  
  findex <- rep(0,p)
  
  findex[xname %in% isFactor] <- 1
  
  designTable <- list('table' = rbind( round(VIF,2),findex[-wi],round(corx,2)) )
  rownames(designTable$table) <- c('VIF','factor',xname[-wi])
  
  designTable$table <- designTable$table[-3,]
  
  if(p == rankx)designTable$rank <- paste('full rank:',rankx,'= ncol(x)')
  if(p < rankx) designTable$rank <- paste('not full rank:',rankx,'< ncol(x)')

  list(VIF = round(VIF,2), correlation = round(corx,2), rank = rankx, p = p,
       isFactor = isFactor, designTable = designTable)
}

.getYscore <- function(yscore, nscore=length(yscore), PLOT=F, cex=1){
  
  ord    <- order(yscore,decreasing=T)
  yo     <- yscore[ord[1:nscore]]
  onames <- names(yo)
  
  if(PLOT){
    par(mfrow=c(1,1),bty='n')
    ylim    <- range(yo)
    ylim[2] <- ylim[1] + 1.2*diff(ylim) 
    xlim    <- c(0,nscore*1.2)
    plot(1:nscore,yo,type='s',xlim=xlim,ylim=ylim)
    text(c(1:nscore),yo,onames,srt=75,pos=4, cex=cex)
  }
  list(scores = yo, index = ord)
}


.fitText2Fig <- function(xx, width=T, fraction=1, cex.max=1){
  
  # returns cex to fit xx within fraction of the current plotting device
  # width - horizontal labels stacked vertically
  #!width - vertical labels plotted horizontally
  
  px <- par('pin')[1]
  py <- par('pin')[2]
  cl <- max( strwidth(xx, units='inches') )
  ch <- strheight(xx, units='inches')[1]*length(xx)  # ht of stacked vector
  
  if(width){              #horizontal labels stacked vertically
    xf <- fraction*px/cl
    yf <- fraction*py/ch
  } else {                #vertical labels plotted horizontally
    xf <- fraction*px/ch
    yf <- fraction*py/cl
  }
  
  cexx <- min(c(xf,yf))
  if(cexx > cex.max)cexx <- cex.max
  cexx
}


.cov2Dist <- function(sigma){ #distance induced by covariance
  
  n <- nrow(sigma)
  matrix(diag(sigma),n,n) + matrix(diag(sigma),n,n,byrow=T) - 2*sigma
}

.distanceMatrix <- function(mat, DIST=F){
  
  # mat is n by m matrix
  # if DIST returns a m by m distance matrix, otherwise corr matrix
  
  if(isSymmetric(mat)){
    if( all(diag(mat) == 1) ){   #now a correlation matrix
      mmm1 <- mat
      if(DIST)mmm1 <- .cov2Dist(mat)
    } else {                      #now a covariance 
      if(DIST){
        mmm1 <- .cov2Dist( mat )
      } else {
        mmm1 <- cor(mat)
      }
    }
  } else  {     # not symmetric
    if(DIST){
      mmm1 <- .cov2Dist( cov(mat) )
    } else {
      mmm1 <- cor(mat)
    }
  }
  mmm1
}



.clusterWithGrid <- function(mat1, mat2=NULL, DIST=F, expand=1,
                             mainLeft=' ', main1 = ' ', main2 = ' ',
                             leftClus=F, rightClus=F, topClus1=F, topClus2=F,
                             leftLab=F, rightLab=F, topLab1=F, topLab2=F,
                             colOrder1=NULL, colOrder2=NULL, rowOrder=NULL, 
                             colCode1 = NULL, colCode2 = NULL, rowCode=NULL,
                             lower1 = F, diag1 = F, lower2 = F, diag2 = F,
                             slim1=NULL, slim2=NULL,
                             horiz1 = NULL,  horiz2 = NULL,
                             vert1 = NULL, vert2 = NULL){
  
  #   layout: mat1 on left, mat2 (if given) on right
  # clusters: left & top or right & top
  #   expand: width of mat1 relative to mat2
  # if cluster analysis is used to order, then mat1 must be symmetric matrix
  # if DIST use distance, otherwise correlation
  
  doneLeft <- done1 <- done2 <- F
  
  nr  <- nrow(mat1)
  nc1 <- ncol(mat1)
  nc2 <- 0
  
  twoMat <- F
  if(!is.null(mat2)){
    twoMat <- T
    nc2 <- ncol(mat2)
 #   mat2 <-  apply(mat2,2,rev) 
    if(nrow(mat2) != nr)stop('matrices must have same no. rows')
  }
  cwide  <- .15
  mg     <- .08
  lg     <- rg <- tg <- mg
  gg <- .24
  
  if(leftLab) lg <- gg
  if(topLab1 | topLab2)  tg <- gg
  if(rightLab)rg <- gg
  
  xwide <- mg
  if(leftLab) xwide <- c(xwide,lg)
  if(leftClus)xwide <- c(xwide,cwide)
  
  xg <- .8
  if(twoMat){
    xg <- expand*nc1/(expand*nc1 + nc2)
    xg <- c(xg,1 - xg)
  }
  xwide <- c(xwide,xg)
  
  if(rightClus)xwide <- c(xwide,cwide)
  if(rightLab) xwide <- c(xwide,rg)
  xwide <- c(xwide,mg)
  xloc <- cumsum(xwide)/sum(xwide)
  
  ywide <- c(mg,.8)
  if(topClus1 | topClus2)ywide <- c(ywide,cwide)
  if(topLab1 | topLab2) ywide <- c(ywide,tg)
  ywide <- c(ywide,mg)
  yloc  <- cumsum(ywide)/sum(ywide)
  
  mmm1 <- .distanceMatrix(mat1, DIST)
  rrr  <- .distanceMatrix(t(mat1), DIST)
  
  dcorr <- distr <- dcor1 <- dist1 <- dcor2 <- dist2 <- NULL
  
  if( DIST){
    distr <- rrr
    dist1 <- mmm1
    labCorr <- rownames(distr)
    labCor1 <- rownames(dist1)
  } else {
    dcorr <- rrr
    dcor1 <- mmm1
    labCorr <- rownames(dcorr)
    labCor1 <- rownames(dcor1)
  }
  
  if(is.null(rowOrder)){    # row clustering
    tmp <- .clusterPlot( dcor=dcorr, dist=distr, PLOT=F)
    rowOrder <- tmp$corder  #order bottom to top
  }
  if(is.null(colOrder1)){      # column clustering
    if(isSymmetric(mat1)){
      colOrder1 <- rowOrder
    } else {
      tmp <- .clusterPlot( dcor=dcor1, dist=dist1, PLOT=F)
      colOrder1 <- tmp$corder
    }
  }
  if(!is.null(mat2) & is.null(colOrder2)){
    mmm2 <- .distanceMatrix(mat2, DIST)
    if( DIST){
      dist2 <- mmm2
    } else {
      dcor2 <- mmm2
    }
    labCor2 <- colnames(mmm2)
    
    tmp <- .clusterPlot( dcor=dcor2, dist=dist2,PLOT=F)
    colOrder2 <- tmp$corder
  }
  
 # mat1 <-  apply(mat1,2,rev)    # organize down to up
  
  if(is.null(rowCode)) rowCode  <- rep('black',nr)
  if(is.null(colCode1))colCode1 <- rep('black',nc1)
  if(is.null(colCode2))colCode2 <- rep('black',nc2)
  
  #######################
  NEW <- add <- F
  xi <- 0:1
  yi <- 1:2
  
  # r1 <- rrr[rowOrder,colOrder1]
  m1 <- mat1[rev(rowOrder),colOrder1]
  if(!is.null(mat2))m2 <- mat2[rev(rowOrder),colOrder2]
  
  ##### lab panel
  if(leftLab){
    xi <- xi + 1
    par(plt=c(xloc[xi],yloc[yi]),bty='n', new=NEW)
    
    plot(c(0,0),c(0,0),col='white',xlim=c(0,1),ylim=c(0,nr),
         xaxt='n',yaxt='n',xlab='',ylab='')
    xl  <- rep(1,nr)
    yl  <- c(1:nr)*nr/diff(par('usr')[3:4])
    cex <- .fitText2Fig(rownames(rrr),fraction=.96)
    text( xl,yl,labCorr[rowOrder],pos=2,cex=cex, 
          col = rowCode[rowOrder])
    NEW <- add <- T
    mtext(mainLeft,2)
    doneLeft <- T
  }
  
  #### cluster panel
  if(leftClus){
    xi <- xi + 1
    par(plt=c(xloc[xi],yloc[yi]),bty='n',  new=NEW)
    .clusterPlot( dcor=dcorr, dist=distr ,main=' ',cex=.2, 
                  colCode=rowCode,
                  LABELS = F,horiz=T, noaxis=T)
    NEW <- add <- T
    if(!doneLeft)mtext(mainLeft,2)
    doneLeft <- T
  }
  
  ######## first grid plot
  
  xi <- xi + 1
  yz <- yi
  
  if(topClus1){
    yz <- yz + 1
    par(plt=c(xloc[xi],yloc[yz]),bty='n',new=NEW)
  
    .clusterPlot( dcor=dcor1, dist=dist1 ,main=' ', 
                  colCode=colCode1[colOrder1],      
                  LABELS = F, horiz=F, noaxis=T, add=T)
    NEW <- add <- T
    if(!topLab1){
      mtext(main1,3)
      done1 <- T
    }
  }
  
  par(plt=c(xloc[xi],yloc[yi]), bty='n', new=NEW)
  if(is.null(slim1))slim1 = quantile(m1,c(.01,.99))
  slim1  <- signif(slim1,1)
  
  tmp    <- .colorSequence(slim1)
  scale  <- tmp$scale
  colseq <- tmp$colseq
  
  ww    <- as.matrix(expand.grid(c(1:nr),c(1:nc1)))  # note reverse order
  mt    <- m1
  
  mask <- lower.tri(mt,diag=!diag1)
  mask <- apply(mask,2,rev)
  
  if(lower1)mt[mask] <- 0
  
  icol <- findInterval(mt[ww],scale,all.inside=T)
  coli <- colseq[icol]
  
  xlim=c(range(ww[,2])); xlim[2] <- xlim[2] + 1
  ylim=c(range(ww[,1])); ylim[2] <- ylim[2] + 1
  
  sides <- cbind( rep(1,nrow(ww)), rep(1,nrow(ww)) )
  plot(0,0,cex=.1,xlab=' ',ylab=' ', col='white',
       xaxt='n',yaxt='n', xlim=xlim, ylim=ylim)
  
  symbols(ww[,2] + .5,nr - ww[,1] + 1 + .5,rectangles=sides,
          fg=coli,bg=coli,inches=F, xlab=' ',ylab=' ',
          xaxt='n',yaxt='n', add=T)
  
  if(!is.null(horiz1)){
    cut <- which(diff(horiz1[colOrder1]) != 0) + 1
    ncc <- length(cut)
    for(i in 1:ncc){
      lines(c(0,cut[i]-2),cut[c(i,i)],lty=2)
    }
    text(rep(1,ncc),cut,2:(ncc+1),pos=3)
  }
  if(!is.null(vert1)){
    cut <- which(diff(vert1[colOrder1]) != 0) + .5
    ncc <- length(cut)
    for(i in 1:ncc){
      lines(cut[c(i,i)],c(cut[i]+2,nc1),lty=2)
    }
    text(cut,rep(nc1,ncc),2:(ncc+1),pos=4)
  }
  
  NEW <- add <- T
  if(!doneLeft)mtext(mainLeft,2)
  doneLeft <- T
  
    
  if(topLab1){
    yz <- yz + 1
    par(plt=c(xloc[xi], yloc[yz]),bty='n', new=NEW)
    plot(c(0,0),c(0,0),col='white',xlim=c(1,nc1) ,ylim=c(0,1),
         xaxt='n',yaxt='n',xlab='',ylab='')
    yl <- rep(0,nc1)
    xl <- .99*c(1:nc1)*(nc1-1)/diff(par('usr')[1:2])
 #   xl <- .95*c(1:nc1)*nc1/diff(par('usr')[1:2])
    cex <- .fitText2Fig(colnames(m1), width=F, fraction=.95)
    text( xl - .5,yl,colnames(m1),pos=4,cex=cex,srt=90,
          col=colCode1[colOrder1])
  }
  if(!done1)mtext(main1,3)
  
  ######## 2nd grid plot
  
  if(twoMat){
    
    xi <- xi + 1
    yz <- yi
 
    if(topClus2){
      yz <- yz + 1
      par(plt=c(xloc[xi],yloc[yz]),bty='n',new=NEW)
      tmp <- .clusterPlot( dcor=dcor2, dist=dist2 ,main=' ', LABELS = F,
                           colCode=colCode2[colOrder2], horiz=F, 
                           noaxis=T, add=T)
      if(!topLab2){
        mtext(main2,3)
        done2 <- T
      }
    }
    
    par(plt=c(xloc[xi],yloc[yi]), bty='n', new=T)
    if(is.null(slim2))slim2 = quantile(m2,c(.01,.99))
    slim2  <- signif(slim2,1)
    
    tmp <- .colorSequence(slim2)
    scale  <- tmp$scale
    colseq <- tmp$colseq
    
    ww    <- as.matrix(expand.grid(c(1:nr),c(1:nc2)))  # note reverse order
    mt    <- m2
    
    if(lower2){
      mask <- lower.tri(mt,diag=!diag1)
      mask <- apply(mask,2,rev)
      mt[mask] <- 0
    }
    
    icol <- findInterval(mt[ww],scale,all.inside=T)
    coli <- colseq[icol]
    
    xlim=c(range(ww[,2])); xlim[2] <- xlim[2] + 1
    ylim=c(range(ww[,1])); ylim[2] <- ylim[2] + 1
    
    sides <- cbind( rep(1,nrow(ww)), rep(1,nrow(ww)) )
    plot(0,0,cex=.1,xlab=' ',ylab=' ',
         col='white',xaxt='n',yaxt='n', xlim=xlim, ylim=ylim)
    
    symbols(ww[,2] + .5,nr - ww[,1] + 1 + .5, rectangles=sides,
            fg=coli, bg=coli, inches=F, xlab=' ',ylab=' ',
            xaxt='n', yaxt='n', add=T)
    
    if(!is.null(horiz2)){
      cut <- which(diff(horiz2[colOrder1]) != 0) + 1
      ncc <- length(cut)
      for(i in 1:ncc){
        xmm <- c(0,cut[i]-2)
        if(!lower2)xmm[2] <- nc2 + 1
        lines(xmm,cut[c(i,i)],lty=2)
      }
      if(lower2) text(rep(1,ncc),cut,2:(ncc+1),pos=3)
      if(!lower2)text(rep(nc2+1,ncc),cut,2:(ncc+1),pos=3)
    }
    if(!is.null(vert2)){
      cut <- which(diff(vert2[colOrder1]) != 0) + .5
      ncc <- length(cut)
      for(i in 1:ncc){
        lines(cut[c(i,i)],c(cut[i]+2,nc1),lty=2)
      }
      text(cut,rep(nc1,ncc),2:(ncc+1),pos=4)
    }
    
    if(topLab2){
      yz <- yz + 1
      par(plt=c(xloc[xi],yloc[yz]),bty='n', new=NEW)
      plot(c(0,0),c(0,0),col='white',xlim=c(1,nc2),ylim=c(0,1),
           xaxt='n',yaxt='n',xlab='',ylab='')
      yl <- rep(0,nc2)
      xl <- .99*c(1:nc2)*(nc2-1)/diff(par('usr')[1:2])
  #    xl <- .95*c(1:nc2)*nc2/diff(par('usr')[1:2])
      cex <- .fitText2Fig(colnames(m2),width=F, fraction=.95)
      text( xl - .2,yl,colnames(m2),pos=4,cex=cex,srt=90, 
            col=colCode2[colOrder2])
    }
    if(!done2)mtext(main2,3)
  }
  
  if(rightClus){
    xi <- xi + 1
    par(plt=c(xloc[xi], yloc[yi]), bty='n', mgp=c(3,1,0), new=NEW)
    tmp <- .clusterPlot( dcor = rrr ,main=' ',cex=.2, REV=T,
                         LABELS = F,horiz=T, noaxis=T)
  }
  
  if(rightLab){
    xi <- xi + 1
    par(plt=c(xloc[xi],yloc[yi]),bty='n', new=NEW)
    plot(c(0,0),c(0,0),col='white',xlim=range(c(0,1)),ylim=c(0,nr),
         xaxt='n',yaxt='n',xlab='',ylab='')
    xl <- rep(0,nr)
    yl <- c(1:nr)*nr/diff(par('usr')[3:4])
    cex <- .fitText2Fig(rownames(m1),fraction=.8)
    text( xl,yl,labCorr[rowOrder],pos=4,cex=cex,
          col=rev(rowCode[rowOrder]))
  }
}


.clusterPlot <- function(dcor=NULL,dist=NULL,main=' ',xlab='Species',
                         method='complete',
                        cex=1, ncluster=2, add=F, REV=F,
                        xlim=NULL, colCode = NULL, horiz=T, textSize=1,
                        LABELS=T, noaxis=F,PLOT=T){  
  
  #dcor is a correlation matrix
  #dist is a distance matrix
  
  if(!is.null(dist)){
    if(!LABELS) rownames(dist) <- colnames(dist) <- NULL
    nn <- nrow(dist)
    diss <- as.dist( dist )
    nr   <- nrow(dist)
  }
  if(is.null(dist)){
    if(!LABELS) rownames(dcor) <- colnames(dcor) <- NULL
    nn <- nrow(dcor)
    diss <- as.dist(.cov2Dist(dcor))
    nr   <- nrow(dcor)
  }
  
  tmp    <- hclust(diss,method)
  corder <- tmp$order
  names(corder) <- tmp$labels
  ctmp   <- cutree(tmp,k=1:ncluster)
 # if(!LABELS)attr(tmp,'labels') <- rep('',nr) 
  
  wclus <- ctmp[,ncluster]
  clusterCol <- NULL
  
  clusterIndex <- ctmp[,ncluster]
  clusterList <- character(0)
  
#  mycols <- mapColors(ncluster)
  
  notLab <- F
  if(is.null(colCode)){
    colF   <- colorRampPalette(c('black','blue','orange','brown','red'))
    mycols <- colF(ncluster)
    notLab <- T
    colCode <- mycols[ctmp[,ncluster]]
    names(colCode) <- rownames(ctmp)
  }
  col.lab <- colCode
  if(!LABELS)col.lab <- rep('white',length(colCode))
  
  colLab <- function(n) {
 
    if(is.leaf(n)) {
      a <- attributes(n)
      attr(n, "nodePar") <- c(a$nodePar, 
                              list(col = col.lab[n[1]],
                                   lab.col = col.lab[n[1]]))
    }
    n
  }
  tdendro <- as.dendrogram(tmp)
  dL      <- dendrapply(tdendro,colLab)
  
  nodePar <- list(cex = .1, lab.cex=textSize)
    leafLab         <- "textlike"
    nodePar$leaflab <-  leafLab
    
  if(!PLOT){
    return(  invisible(list( clusterList = clusterList, colCode = colCode, 
                    clusterIndex = clusterIndex,
                    corder = corder) ) )
  }
    
    if(horiz){
      if(is.null(xlim))xlim <- c(attr(dL,'height'),0)
      if(REV)xlim <- rev(xlim)
    }
    
  axes <- T
  if(noaxis)axes <- F

  new <- F
  if(add)new <- T
 # par(new = new,cex=textSize)
  tmp <- plot( dL,nodePar=nodePar, horiz=horiz, xlim=xlim, 
               axes = axes)
  if(!LABELS & !notLab){
    
    col <- colCode[corder]
    pvar <- par('usr')
    
    wi <- abs(diff(pvar[1:2])/10)
    hi <- abs(diff(pvar[3:4])/10)
    
    if(horiz){
      xx <- rep(pvar[2],nn)
      yy <- 1:nn
      rec <- cbind( rep(wi,nn), rep(1,nn) )
      symbols(xx,yy,rectangles=rec,fg=col, bg=col, inches=F, add=T)
    } else {
      xx <- 1:nn
      yy <- rep(pvar[3],nn)
      rec <- cbind( rep(1,nn), rep(hi,nn) )
      symbols(xx,yy,rectangles=rec,fg=col, bg=col, inches=F, add=T)
    }
  }
  
  title(main)
  
  invisible(list( clusterList = clusterList, colCode = colCode, 
                  clusterIndex = clusterIndex,
                  corder = corder) )

}
.colorLegend <- function(xx,yy,ytick=NULL,
                         scale=seq(yy[1],yy[2],length=length(cols)),
                        cols,labside='right', text.col=NULL,
                        bg=NULL,endLabels=NULL){  
  # xx = (x1,x2), y = (y1,y2)
  # bg is color of border

  nn <- length(scale)
  ys <- seq(yy[1],yy[2],length=nn)

  for(j in 1:(length(scale)-1)){
    rect(xx[1],ys[j],xx[2],ys[j+1],col=cols[j],border=NA)
  }
  if(!is.null(bg))rect(xx[1],yy[1],xx[2],yy[2],border=bg,lwd=3)
  if(!is.null(ytick)){
    
    ys <- diff(yy)/diff(range(ytick))*ytick
    yt <- ys - min(ys) + yy[1]
    
    for(j in 1:length(yt)){
      lines(xx,yt[c(j,j)])
    }
  }
  if(!is.null(endLabels)){ 
    cx <- cols[c(1,nn)]
    if(!is.null(text.col))cx <- text.col
    if(!is.null(text.col))cx <- text.col
    if(labside == 'right')text(diff(xx)+c(xx[2],xx[2]),yy,endLabels,col=cx)
    if(labside == 'left')text(c(xx[1],xx[1]),yy,endLabels,pos=2,col=cx)
  }
}
.corMatCI <- function(rmat, nmat, alpha = .05){
  
  # rmat   - m by m correlation matrix
  # nmat   - m by m matrix of counts
  # lo, hi - lower, upper CI
  # sig    - rmat outside CI
  
  m <- nrow(rmat)
  if(!is.matrix(nmat))nmat <- matrix(nmat,m,m)
  
  lo <- hi <- sr <- nmat*0
  ii <- lower.tri(nmat)
  
  fz <- .5*log( (1 + rmat[ii])/(1 - rmat[ii]) )
  se <- 1/sqrt(nmat[ii] - 3)
  ds <- se*qnorm(1 - alpha/2)
  lo[ii] <- .fishz2r(fz - ds)
  hi[ii] <- .fishz2r(fz + ds)
  sr[ (rmat < lo | rmat > hi) & ii] <- 1
  
  list(lo = lo, hi = hi, sig = sr)
}

.capFirstLetter <- function(xx) {     
  
  #capiltalize first letter of every word
  
  s <- unlist(strsplit(xx, " "))
  s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
             sep = "", collapse = " ")
  unlist(strsplit(s, " "))
}

.lowerFirstLetter <- function(xx){
  s <- unlist(strsplit(xx, " "))
  s <- paste(tolower(substring(s, 1, 1)), substring(s, 2),
             sep = "", collapse = " ")
  unlist(strsplit(s, " "))
}



.colorSequence <- function(slim, colorGrad=NULL, ncol=200){  
  
  # generate color sequence with white for zero
  # slim is scale from min to max
  # used in .corPlot
  
  if(is.null(colorGrad)){
    colorSeq <- c('darkblue','darkblue','blue',
                  'green','white',
                  'yellow','red','brown','brown')
    colorGrad   <- colorRampPalette(colorSeq)
  }
  
  colseq <- colorGrad(ncol)
  
  if(slim[1] < 0 & slim[2] > 0){  #put white at zero
    dp <- slim[2] - 0
    dm <- 0 - slim[1]
    ncol <- 200
    
    colseq <- colorGrad(ncol)
    if(dp < dm)colseq <- colseq[101 + c(-100:round(dp/dm*100))]
    if(dp > dm)colseq <- colseq[ round((1 - dm/dp)*100):200 ]
    ncol  <- length(colseq)
  }
  scale <- seq(slim[1],slim[2],length.out=ncol)
  return( list(colseq = colseq, scale = scale ) )
}

.corPlot <- function(cmat,slim=NULL,PDIAG=F,plotScale=1,
                    makeColor=NULL,textSize=NULL,
                    textCol = rep('black',nrow(cmat)), 
                    corLines=T,tri='lower',colorGrad = NULL,
                    cex=1, specLabs = T, squarePlot = T,LEGEND = T,
                    widex=5.5,widey=6.5,add=F,new=T){  #correlation or covariance matrix

  # makeColor - list of matrices of indices for boxes
  #   names of matrices are colors
  # if(PDIAG)diag(cmat) <- 0
  # tri - 'lower','upper', or 'both'
  # colorGrad - constructed with colorRampPalette()
  # squarePlot makes symbols square
  # new means NOT NEW 

  if(is.null(slim))slim = quantile(cmat,c(.01,.99))
  slim  <- signif(slim,1)
  
  if(tri == 'upper')cmat[lower.tri(cmat)] <- 0
  if(tri == 'lower')cmat[upper.tri(cmat)] <- 0

  dy  <- nrow(cmat)
  dx  <- ncol(cmat)
  d <- dx
  xtext <- rep(c(1,100),dx/2)
  if(length(xtext) < d)xtext <- c(xtext,1)

  if(d < 20)xtext <- xtext*0 + 1

  xtext <- xtext*0 + 1
  
  if(!is.null(colorGrad)){
    ncol  <- 200
    colseq <- colorGrad(ncol)
    scale  <- seq(slim[1],slim[2],length.out=ncol)
  } else {
    tmp <- .colorSequence(slim, colorGrad)
    scale  <- tmp$scale
    colseq <- tmp$colseq
  }
  
  ww   <- as.matrix(expand.grid(c(1:dy),c(1:dx)))  # note reverse order

  if(tri == 'upper'){
    ww  <- ww[ww[,1] <= ww[,2],]
    ww  <- ww[order(ww[,1]),]
  }
  if(tri == 'lower'){
    ww  <- ww[ww[,1] >= ww[,2],]
    ww  <- ww[order(ww[,1]),]
  }

  icol <- findInterval(cmat[ww],scale,all.inside=T)
  coli <- colseq[icol]

  if(PDIAG)coli[ww[,1] == ww[,2]] <- 'white'
  
  ss <- max(c(dx,dy))/5/plotScale
  
  if(squarePlot).mapSetup(c(0,dx),c(0,dy),scale=ss,
                          widex=widex,widey=widey)
  
  if(squarePlot){
    symbols(ww[,2],dy - ww[,1] + 1,squares=rep(1,nrow(ww)),
            xlim=c(0,dx+4),ylim=c(0,dy+4),
            fg=coli,bg=coli,inches=F,xlab=' ',ylab=' ',xaxt='n',yaxt='n',
            add=add)
  } else {
    sides <- cbind( rep(1,nrow(ww)), rep(1,nrow(ww)) )
    symbols(ww[,2],dy - ww[,1] + 1,rectangles=sides,
            xlim=c(0,dx+4),ylim=c(0,dy+4),
            fg=coli,bg=coli,inches=F,xlab=' ',ylab=' ',xaxt='n',yaxt='n',
            add=add)
  }
  
  if(!is.null(makeColor)){
    
    for(k in 1:length(makeColor)){
      mm <- makeColor[[k]]
      if(length(mm) == 0)next
      if(tri == 'upper')mm <- mm[mm[,1] <= mm[,2],]
      if(tri == 'lower')mm <- mm[mm[,1] >= mm[,2],]
      ss <- matrix(0,dy,dx)
      ss[mm] <- 1
      wk <- which(ss[ww] == 1)
      ccc <- names(makeColor)[[k]]
      symbols(ww[wk,2],dy - ww[wk,1]+1,squares=rep(1,length(wk)),
              fg=ccc,bg=ccc,inches=F,xaxt='n',yaxt='n',add=T)
    }
  }
  
  ncolor <- length(unique(textCol))
  
  ll <- 1/d + 1
  
  if(tri == 'lower'){
    for(kk in 1:d){
      kb <- kk - .5
      ke <- d - kk + .5
      
      if(corLines){
        if(kk <= d)lines(c(kb,kb),c(0,ke),col='grey',lwd=1.5) #vert
        if(kk > 1){
          lines( c(.5,kb),c(ke,ke),col='grey',lwd=1.5)        #horizontal
          lines(c(kb,kb+.5),c(ke,ke+.5),col='grey',lwd=1.5)   #diagonal
        }
      }
      if(!specLabs & ncolor > 1){
        xp <- c(kb, kb, kb + ll + .5, kb + ll + 1.5, kb + 1)
        yp <- c(ke, ke + 1, ke + ll + 1.5, ke + ll + .5, ke)
        polygon(xp, yp, border = textCol[kk], col = textCol[kk])
      }
    }
  }
  rect(0,-1,d+1,.5,col='white',border=NA)
  
  if(is.null(textSize))textSize <- exp(-.02*ncol(cmat))
  labels   <- rev(rownames(cmat))
  if(!specLabs)labels <- F
  
  if(tri == 'lower' & specLabs)text( c(d:1)+.1*xtext, c(1:d)+.5, 
                                     rev(colnames(cmat)),pos=4,srt=45,
                                     col = rev(textCol), cex=textSize)
  
  if(tri == 'both'){
    labels   <- rev(rownames(cmat))
    par(las = 1)
 #   tmp <- axis(side=2,at=c(1:dy),labels=labels,tick=F,lwd=0,pos=.5,
 #        cex.axis=textSize, col = rev(textCol))
    
    .yaxisHorizLabs( labels, at=c(1:length(labels)), xshift=.05,
                                 col = textCol, pos=2)
    
    par(las = 0)
    
    if(specLabs){
      text( c(dx:1)-.1*xtext, xtext*0+dy+.8, rev(colnames(cmat)),
                      pos=4, srt=55, col = rev(textCol), cex=textSize)
    } else {
      sides <- cbind( rep(1,dx),rep(1/dy,dx) )
      symbols(1:dx,rep(1+dy,dx),rectangles=sides,
              fg=textCol,bg=textCol,
              add=T)
    } 
      
  }
  
  labside <- 'left'
  
  wk <- which(scale >= slim[1] & scale <= slim[2]) 
  
  px <- par('usr')
  xs <- .01*diff(px[1:2])
  midx <- .95*mean( c(dx,px[2]) )
  
  yx <- c(.2*dy,.2*dy + .35*dy)
  
  if(LEGEND).colorLegend(c(midx-xs,midx+xs),yx,ytick=c(slim[1],0,slim[2]),
                        scale[wk],cols=colseq[wk],labside=labside,
                        endLabels=range(slim),text.col='black')
}
.cor2Cov <- function(sigvec,cormat){ #correlation matrix and variance vector to covariance
  
  d <- length(sigvec)
  s <- matrix(sigvec,d,d)
  cormat*sqrt(s*t(s))
}

.cov2Cor <- function(covmat, covInv = NULL){  
  # covariance matrix to correlation matrix
  # if covInv provided, return inverse correlation matrix

  d    <- nrow(covmat)
  di   <- diag(covmat)
  s    <- matrix(di,d,d)
  cc   <- covmat/sqrt(s*t(s))
  
  if(!is.null(covInv)){
    dc <- diag(sqrt(di))
    ci <- dc%*%covInv%*%dc
    return(ci)
  }
  cc
}

.cov2Dist <- function(sigma){ #distance induced by covariance
	
	n <- nrow(sigma)
	matrix(diag(sigma),n,n) + matrix(diag(sigma),n,n,byrow=T) - 2*sigma
}
.dMVN <- function(xx,mu,smat=NULL,sinv=NULL,log=F){          #MVN density for mean 0
  
  if(!is.matrix(xx))xx <- matrix(xx,1)
  if(!is.matrix(mu))mu <- matrix(mu,1)
  
  xx <- xx - mu
  
  if(!is.null(sinv)){
    distval <- diag( xx%*%sinv%*%t(xx) )
    ev      <- eigen(sinv, only.values = TRUE)$values
    logd    <- -sum(log(ev))
  }
  
  if(is.null(sinv)){
    testv <- try(chol(smat),T)
    if(inherits(testv,'try-error')){
       tiny  <- min(abs(xx))/100 + 1e-5
       smat  <- smat + diag(diag(smat + tiny))
       testv <- try(chol(smat),T)
    }
    covm    <- chol2inv(testv)
    distval <- rowSums((xx %*% covm) * xx)
    ev      <- eigen(smat, only.values = TRUE)$values 
    logd    <- sum(log( ev ))
  }

    z <- -(ncol(xx) * log(2 * pi) + logd + distval)/2
    if(!log)z <- exp(z)
    z
}

.directIndirectCoeffs <- function( snames, xvector, chains, MEAN = T,
                                   factorList = NULL, keepNames, omitY,
                                   sdScaleY = F, sdScaleX, standX, 
                                   otherpar = NULL, REDUCT = F, ng, burnin,
                                   nsim = 50){
  
  # if MEAN, then use means, otherwise median
  # indirect do not change with x, can choose not to calculate
  #          - a list of vectors, one for each multilevel factor, 
  #            where hostNames appear in colnames of bchain
  #indirFrom - effect from all others
  #indirTo   - effect on all others
  
  if(is.matrix(xvector))stop('xvector must be a row vector with variable names')
  
  xnames <- names(xvector)
  
  N      <- otherpar$N
  r      <- otherpar$r
  bchain <- chains$bgibbs
  schain <- chains$sgibbs
  sigErrGibbs <- kchain <- NULL
  if(REDUCT){
    kchain      <- chains$kgibbs
    sigErrGibbs <- chains$sigErrGibbs
  }
  
  ns <- nsim
  simIndex <- sample(burnin:ng,ns,replace=T)
  
  if(sdScaleY){
    tmp <- .expandSigmaChains(snames, sgibbs = schain, otherpar = otherpar, 
                              simIndex = simIndex, sigErrGibbs, kchain, 
                              REDUCT)
  #  bchain <- tmp$chainList$bchainCor    # standardized
    if(REDUCT)kchain <- kchain[simIndex,]
    schain <- schain[simIndex,]          # not standardized
    sigErrGibbs <- sigErrGibbs[simIndex]
  } else {
    bchain <- bchain[simIndex,]
    schain <- schain[simIndex,]
  }
  
  if(length(factorList) > 0){
    factorNames <- factorList
    for(j in 1:length(factorList)){
      tmp <- matrix( unlist(strsplit(factorList[[j]],names(factorList)[j])),
                     ncol=2,byrow=T)[,2]
      tmp[nchar(tmp) == 0] <- paste(names(factorList)[j],c(1:length(tmp)),
                                    sep='')
      factorNames[[j]] <- tmp
    }
  }
  
  S <- S1 <- length(snames)
  sindex <- c(1:S)
  knames <- snames
  
  nc <- nrow(bchain)
  
  gs <- 1:nrow(bchain)
  
  if(length(omitY) > 0){
    wob <- grep(paste(omitY,collapse="|"),colnames(bchain))
    bchain[,wob] <- 0
    sindex <- sindex[!snames %in% omitY]
    knames <- snames[sindex]
    S1     <- length(knames)
  }
  
  nspec <- length(snames)
 
  ww   <- grep(':',xnames)
  main <- xnames
  if(length(ww) > 0)main <- xnames[-ww]
  main <- main[main != 'intercept']
  int  <- unique( unlist( strsplit(xnames[ww],':') ) ) 
  
  mainEffect <- matrix(NA,nspec,length(main))
  colnames(mainEffect) <- main
  rownames(mainEffect) <- snames
  intEffect  <- dirEffect <- indEffectTo <- mainEffect
  mainSd <- dirSd <- intSd <- indSdTo <- mainEffect 
  
  maxg <- length(main)*length(sindex)*length(gs)
  pbar <- txtProgressBar(min=1,max=maxg,style=1)
  ig   <- 0
  
  for(j in 1:length(main)){
    
    ttt <- .interactionsFromGibbs(mainx=main[j], bchain=bchain,
                                  specs=snames, xmnames=names(xvector), 
                                  xx=xvector, omitY = omitY, sdScaleX=F, 
                                  standX)
    maine   <- ttt$main
    inter   <- ttt$inter  #
    indirTo <- maine*0
    direct  <- maine + inter  # already standardized for X
    
    if(MEAN){
      dmain  <- colMeans(maine)
      inte   <- colMeans(inter)
      dir    <- colMeans(direct)
    } else {
      dmain  <- apply(maine,2,median)
      inte   <- apply(inter,2,median)
      dir    <- apply(direct,2,median)
    }
    
    mainEffect[sindex,j] <- dmain
    intEffect[sindex,j]  <- inte
    dirEffect[sindex,j]  <- dir
    
    mainSd[sindex,j] <- apply(maine,2,sd)
    intSd[sindex,j]  <- apply(inter,2,sd)
    dirSd[sindex,j]  <- apply(direct,2,sd)
    
    for(g in gs){
      
      if(REDUCT){
        Z  <- matrix(schain[g,],N,r)
        ss <- .expandSigma(sigErrGibbs[g], S, Z = Z, kchain[g,], 
                           REDUCT = T)[sindex,sindex]
        if(sdScaleY)cc <- .cov2Cor(ss)
      } else {
        ss <- .expandSigma(schain[g,], S = S, REDUCT = F)[sindex,sindex]
        if(sdScaleY)cc <- .cov2Cor(ss)
      }
      
      for(s in 1:length(sindex)){
        
        if(REDUCT){
          si <- .invWoodburryArma(sigErrGibbs[g], Z[kchain[g,sindex[-s]],])
          if(sdScaleY){
            dc <- diag(sqrt(diag(ss)))[-s,-s]
            ci <- dc%*%si%*%dc
          }
        } else {
          si <- solve(ss[-s,-s])
          if(sdScaleY)ci <- solve(cc[-s,-s])
        }
        
        if(!sdScaleY){
          sonk <- ss[drop=F,s,-s]
          e2   <- sonk%*%si%*%direct[g,-s]
        } else {
          sonk <- cc[drop=F,s,-s]
          e2   <- sonk%*%ci%*%direct[g,-s]      # correlation scale
        }
        indirTo[g,s] <- e2
        
        ig <- ig + 1
        setTxtProgressBar(pbar,ig)
        
      } ##############
    }
    
    if(MEAN){
      indirectTo   <- colMeans(indirTo[gs,])
    } else {
      indirectTo   <- apply(indirTo[gs,],2,median)
    }
    indEffectTo[sindex,j]   <- indirectTo
    indSdTo[sindex,j]       <- apply(indirTo[gs,],2,sd)
  } ######################################
  
  if(!is.null(keepNames)){
    wk <- which(rownames(mainEffect) %in% keepNames)
    mainEffect <- mainEffect[wk,]
    intEffect <- intEffect[wk,]
    dirEffect <- dirEffect[wk,]
    indEffectTo   <- indEffectTo[wk,]
    mainSd    <- mainSd[wk,]
    dirSd     <- dirSd[wk,]
    indSdTo   <- indSdTo[wk,]
  }
  
 # if(!is.null(standardX)){
 #   wx <- match(colnames(mainEffect),names(standardX))
 #   sx <- matrix(standardX[wx],nrow(mainEffect),length(wx),byrow=T)
 #   mainEffect <- mainEffect*sx
 #   intEffect  <- intEffect*sx
 #   dirEffect  <- dirEffect*sx
 #   indEffectTo <- indEffectTo*sx
 # }
  
  list(mainEffect = mainEffect, intEffect = intEffect, dirEffect = dirEffect,
       indEffectTo = indEffectTo, mainSd = mainSd, dirSd = dirSd,
       intSd = intSd, indSdTo = indSdTo)
}

.interactionsFromGibbs <- function(mainx,bchain,specs,xmnames=names(xx),
                                   xx=colMeans(xx), omitY=NULL, sdScaleX, 
                                   standX){
  
  # returns main effects and interactions for variable named main
  # xx are values of covariates to condition on
  # mainx is the name of a main effect
  
  if(length(omitY) > 0){
    wob <- numeric(0)
    for(k in 1:length(omitY)){
      wob <- c(wob, grep(omitY[k],colnames(bchain)))
    }
    bchain[,wob] <- 0
    specs <- specs[!specs %in% omitY]
  }
  
  ww   <- grep(':',xmnames)
  int  <- unique( unlist( strsplit(xmnames[ww],':') ) ) 
  int  <- int[int != mainx]
  
  xj <- paste(mainx,specs,sep='_')
  wj <- which(colnames(bchain) %in%  xj)
  if(length(wj) == 0){
    xj <- paste(specs,mainx,sep='_')
    wj <- which(colnames(bchain) %in%  xj)
  }
  
  maine <- bchain[,xj]
  inter <- maine*0
  
  m1 <- paste(mainx,':',sep='')
  m2 <- paste(':',mainx,sep='')
  i1 <- grep( m1,xmnames )
  i2 <- grep( m2,xmnames )
  
  if(sdScaleX)maine <- maine*standX[mainx,'xsd']  #standardize main effect
  
  if( length(i1) > 0 ){
    
    ww <- match(unlist( strsplit(xmnames[i1],m1) ),xmnames)
    ox <- xmnames[ww[is.finite(ww)]]
    for(kk in 1:length(i1)){
      xi <- paste(xmnames[i1[kk]],specs,sep='_')
      wi <- which(colnames(bchain) %in%  xi)
      if(length(wi) == 0){
        xi <- paste(specs,xmnames[i1[kk]],sep='_')
        wi <- which(colnames(bchain) %in%  xi)
      }
      xik   <- xx[ox[kk]]
      bik   <- bchain[,xi]
      if(sdScaleX){
        xik <- (xik - standX[ox[kk],'xmean'])/standX[ox[kk],'xsd']
        bik <- bik*standX[mainx,'xsd']*standX[ox[kk],'xsd']
      }
      inter <- inter + bik*xik
    }
  }
  
  if( length(i2) > 0 ){
    
    ww <- match(unlist( strsplit(xmnames[i2],m2) ),xmnames)
    ox <- xmnames[ww[is.finite(ww)]]
    for(kk in 1:length(i2)){
      xi <- paste(xmnames[i2[kk]],specs,sep='_')
      wi <- which(colnames(bchain) %in%  xi)
      if(length(wi) == 0){
        xi    <- paste(specs,xmnames[i2[kk]],sep='_')
        wi <- which(colnames(bchain) %in%  xi)
      }
      xik   <- xx[ox[kk]]
      bik   <- bchain[,xi]
      if(sdScaleX){
        xik <- (xik - standX[ox[kk],'xmean'])/standX[ox[kk],'xsd']
        bik <- bik*standX[mainx,'xsd']*standX[ox[kk],'xsd']
      }
      inter <- inter + bik*xik
    }
  }
  list(main = maine, inter = inter)
}


.stackedBoxPlot <- function( stackList, stackSd=character(0),
                            ylim=NULL,sortBy = NULL, barnames=NULL,
                            col=rep(NULL,length(stackList)),
                            border=rep(NA,length(stackList)),
                            decreasing=T, nsd=1.96, cex=1,
                            legend=NULL, scaleLegend=.1){
  
  # sortBy - if length 1 indicates which variable in stackList to sort by
  #        - if a vector it is the order to plot
  # nds    - no. standard deviations for whiskers
  
  nn  <- length(stackList)
  ord <- c(1:length(stackList[[1]]))
  nx  <- length(ord)
  
  xx <- 0:(nx-1)
  
  if(is.null(ylim)){
    
    ymax <- ymin <- 0
    
    for(j in 1:nn){
      ymax <- ymax + max( c(0,stackList[[j]]),na.rm=T )
      ymin <- ymin + min( c(0,stackList[[j]]),na.rm=T )
    }
    
    ylim <- c(ymin,ymax)
    
    yscale <- diff(ylim,na.rm=T)*.4
    ylim[1] <- ylim[1] - yscale
    ylim[2] <- ylim[2] + yscale
  }
  
  if(!is.null(sortBy)){
    
    if(length(sortBy) > 1){
      ord <- sortBy
    } else {
      ord <- order( stackList[[sortBy]], decreasing = decreasing)
    }
    if(!is.null(barnames))barnames <- barnames[ord]
  }
  
  dy   <- diff(ylim)
  xlim <- c(0,1.2*length(ord))

  
  add <- F
  
  offset <- offsetPos <- offsetNeg <- rep(0,length(stackList[[1]]))
  
  if(is.null(col))col <- c(1:nn)
  
  for(j in 1:nn){
    
    xj <- stackList[[j]][ord]
    names(xj) <- NULL
    
    wp <- which(xj > 0)     # increase
    wn <- which(xj < 0)     # decrease
    
    offset[wp] <- offsetPos[wp]
    offset[wn] <- offsetNeg[wn]
    
    hj <- xj 
    
    barplot(height= hj,offset=offset,xlim=xlim,ylim=ylim,
            col=col[j],border=border[j],add=add)
    
    ww <- grep(names(stackList)[j],names(stackSd))
    if(length(ww) > 0){
      xz <- xx + .5
      xz <- xz*1.2
      
      tall <-  nsd*stackSd[[ww]]
      y1   <-  hj + offset + tall
      y2   <-  hj + offset - tall
      
      for(i in 1:length(ord)){
        lines(xz[c(i,i)],c(y1[i],y2[i]),lwd=6,col='white')
        lines(c(xz[i]-.1,xz[i]+.1),y1[c(i,i)],lwd=6,col='white')
        lines(c(xz[i]-.1,xz[i]+.1),y2[c(i,i)],lwd=6,col='white')
        
        lines(xz[c(i,i)],c(y1[i],y2[i]),lwd=2,col=col[j])
        lines(c(xz[i]-.1,xz[i]+.1),y1[c(i,i)],lwd=2,col=col[j])
        lines(c(xz[i]-.1,xz[i]+.1),y2[c(i,i)],lwd=2,col=col[j])
      }
    }
    
    if(j == 1)add <- T
    
    offsetPos[wp] <- offsetPos[wp] + hj[wp]
    offsetNeg[wn] <- offsetNeg[wn] + hj[wn]
    
    if(j == nn & !is.null(barnames)){
      
      xall <- par('usr')[1:2]
      xtic <- c(1:nx)*(diff(xall) - 1)/nx - .8
      
      yy <- offsetPos + .2*dy
      pos <- yy*0 + 1
      wl <- which(abs(offsetNeg) < abs(offsetPos))
      yy[wl] <- offsetNeg[wl] - .2*dy
      pos[wl] <- 4
      text(xtic,yy,barnames,srt=90.,pos=pos,cex=cex)
    }
  } 
  
  if(!is.null(legend)){
    
    dy <- diff(ylim)*scaleLegend
    dx <- 1.2
    x1 <- length(ord)*.02 + 1
    y1 <- ylim[1]
    pos <- 4
    if(legend == 'topright'){
      x1  <- length(ord)
      y1  <- ylim[2]
      dy  <- -dy
      dx <- -dx
      pos <- 2
    }
    if(legend == 'topleft'){
      y1  <- ylim[2]
      dy  <- -dy
    }
    if(legend == 'bottomright'){
      x1  <- length(ord)
      dx <- -dx
      pos <- 2
    }
    for(j in 1:length(stackList)){
      y2 <- y1 + dy
      rect(x1,y1,x1 + 1,y2,col=col[j],border=border[j])
      y1 <- y2
      colj <- col[j]
      if(colj == 'white')colj <- border[j]
      text(x1 + dx,y1 - dy/2,names(stackList)[[j]],col=colj,pos=pos,cex=cex)
    }
  }
  
  invisible( ord )
}  



.fishz2r <- function(z){ (exp(2*z) - 1)/(exp(2*z) + 1) }

.getScoreNorm <- function(x,mu,xvar){  #Gneiting/ Raftery proper scoring rule

  #outcome x, prediction mean variance (mu, xvar)

  - ( (x - mu)^2)/xvar - log(xvar)

}

.gjamBaselineHist <- function(y1, bins=NULL, nclass=20){
  
  # add histogram to base of current plot
  
  if(!is.null(bins)){
    hh <- hist(y1,breaks=bins,plot=F)
  } else {
    hh <- hist(y1,nclass=nclass,plot=F)
  }
  
  xvals <- rep(hh$breaks,each=2)
  yvals <- rep(hh$density,each=2)
  
  nb    <- length(hh$breaks)
 # xvals <- c( xvals, hh$breaks[nb] )
  yvals <- c( 0, yvals, 0)
  
  rbind(xvals,yvals)
}

.gjamCensorSetup <- function(y,w,z,plo,phi,wm,censorMat){
  
  nc <- ncol(censorMat)
  br <- numeric(0)
  nk <- length(wm)
  n  <- nrow(y)
  
  zk <- y[,wm]*0
  blast <- -Inf
  
  for(j in 1:nc){
    
    valuej <- censorMat[1,j]
    bj     <- censorMat[2:3,j]
    names(bj) <- paste('c-',names(bj),j,sep='')
    
    if(j > 1){
      if(censorMat[2,j] < censorMat[3,j-1] )
        stop('censor intervals must be unique')
      if(bj[1] == br[length(br)])bj <- bj[2]
    }
    br <- c(br,bj)
    nb <- length(br)
    
    zk[ y[,wm] > blast & y[,wm] < bj[1] ] <- nb - 2
    zk[ y[,wm] == valuej ] <- nb - 1
    blast <- br[length(br)]
  }
  
  if(nc == 1){
    zk[zk == 0] <- 2
    br <- c(br,Inf)
  }
  
  zk[zk == 0] <- 1
  br <- matrix(br,nk,length(br),byrow=T)
  
  censk    <- which(y[,wm] %in% censorMat[1,])
  z[,wm]   <- zk
  
  tmp   <- .gjamGetCuts(z,wm)
  cutLo <- tmp$cutLo
  cutHi <- tmp$cutHi
  
  plo[,wm] <- br[cutLo]
  phi[,wm] <- br[cutHi]
  
  ww <- which(plo[,wm,drop=F] == -Inf,arr.ind=T)
  if(length(ww) > 0){
    mm <- apply(w[,wm],2,max)
    plo[,wm][ww] <- -10*mm[ww[,2]]
  }
  
  tmp <-  .tnorm(nk*n,plo[,wm],phi[,wm],w[,wm],1)
  
  w[,wm][censk] <- tmp[censk]
  
  imat <- w*0                    #location in full matrix
  imat[,wm][censk] <- 1
  censValue <- which(imat == 1)
  
  list(w = w, z = z, cutLo = cutLo, cutHi = cutHi, plo = plo, phi = phi, censValue = censValue,
       breakMat = br)
}

.gjamCuts2theta <- function(tg,ss){   # variance to correlation scale
  
  if(length(ss) == 1)return(tg/sqrt(ss))
  nc   <- ncol(tg)
  sr   <- nrow(ss)
  tg/matrix( sqrt(diag(ss)),sr,nc)
}

.gjamGetCuts <- function(zz,wk){
  
  nk <- length(wk)
  n  <- nrow(zz)

  cutLo <- cbind( rep(1:nk,each=n), as.vector(zz[,wk]) )
  cutHi <- cbind( rep(1:nk,each=n), as.vector(zz[,wk]) + 1 )
  
  list(cutLo = cutLo, cutHi = cutHi)
}
.gjamGetTypes <- function(typeNames=NULL){
  
  TYPES <- c('CON','PA','CA','DA','CAT','FC','CC','OC')
  FULL  <- c('continuous','presenceAbsence','contAbun','discAbun',
             'categorical','fracComp','countComp','ordinal')
  LABS  <- c('Continuous','Presence-absence','Continuous abundance',
             'Discrete abundance',
             'Categorical','Fractional composition','Count composition','Ordinal')
  
  if(is.null(typeNames)){
    names(FULL) <- TYPES
    return( list(typeCols = NULL, TYPES = TYPES, typeFull = FULL, labels = LABS ) )
  }
  
  S        <- length(typeNames)
  typeCols <- match(typeNames,TYPES)
  
  ww <- which(is.na(typeCols))
  if(length(ww) > 0)stop( paste('type code error',typeNames[ww],sep=' ') )
  
  list(typeCols = typeCols, TYPES = TYPES, typeFull = FULL[typeCols],
       typeNames = typeNames, labels = LABS[typeCols])
}

.gjamHoldoutSetup <- function(holdoutIndex,holdoutN,n){
  
  #holdout samples
  if(length(holdoutIndex) > 0)holdoutN <- length(holdoutIndex)
  if(holdoutN > (n/5))stop('too many holdouts')
  
  inSamples <- c(1:n)
  if(holdoutN > 0){
    if(length(holdoutIndex) == 0)holdoutIndex <- sort( sample(n,holdoutN) )
    inSamples <- inSamples[-holdoutIndex]
  }
  nIn <- length(inSamples)
  
  list(holdoutIndex = holdoutIndex, holdoutN = holdoutN, inSamples = inSamples, nIn = nIn)
}
.gjamMissingValues <- function(x,y){
  
  n <- nrow(x)
  
  # missing values in x
  xmiss  <- which(!is.finite(x),arr.ind=T)
  nmiss  <- length(xmiss)
  missX  <- missX2 <- xprior <- yprior <- numeric(0)
  
  xbound <- apply(x,2,range,na.rm=T)
  
  if(nmiss > 0){         #initialize missing values with means
    xmean    <- colMeans(x,na.rm=T)
    x[xmiss] <- xmean[xmiss[,2]]
    xprior   <- x[xmiss]
    nmiss    <- nrow(xmiss)
    fmiss    <- round(100*nmiss/length(x[,-1]),1)
    warning( paste(nmiss,' values (',fmiss,'%) missing in x imputed'), sep='' )
    missX <- missX2 <- rep(0,nmiss)
  }
  
  # rare y
  tmp  <- gjamTrimY(y,minObs=0,OTHER=F)
  wzo  <- which(tmp$nobs == 0)
  if(length(wzo) > 0){
    stop( 'remove from ydata types never present:',
          paste0(names(wzo),collapse=', '))
  }
  fobs <- tmp$nobs/n
  wlo  <- which(fobs < .01)
  if(length(wlo) > 0){
    flo <- paste0(names(fobs)[wlo],collapse=', ')
    warning( paste('present in < 1% of obs:',flo) )
  }
  
  # missing values in y
  ymiss <- which(!is.finite(y),arr.ind=T)
  mmiss <- nrow(ymiss)
  missY <- missY2 <- numeric(0)
  
  if(mmiss > 0){         #initialize missing values with means by TYPEs
    ymean    <- colMeans(y,na.rm=T)
    y[ymiss] <- ymean[ymiss[,2]]
    yprior   <- jitter(y[ymiss])
    fmiss    <- round(100*mmiss/length(y),1)
    warning( paste(mmiss,'values (',fmiss,'%) missing in y are imputed') )
    mmiss <- nrow(ymiss)
    missY <- missY2 <- rep(0,mmiss)
  }
  
  list(xmiss = xmiss, xbound = xbound, missX = missX, missX2 = missX2,
       ymiss = ymiss, missY = missY, xprior = xprior, yprior = yprior)
}
.gjamPlotPars <- function(type='CA',y1,yp,censm=NULL){
  
  if(!is.matrix(y1))y1 <- matrix(y1)
  if(!is.matrix(yp))yp <- matrix(yp)
  
  n       <- nrow(y1)
  nk      <- ncol(y1)
  nbin    <- NULL
  nPerBin <- max( c(10,n*nk/15) )
  breaks  <- NULL
  xlimit  <- quantile(y1,c(.0001,.995))
  ylimit  <- quantile(yp,c(.005,.998),na.rm=T)
  vlines  <- NULL
  wide    <- NULL
  MEDIAN  <- T
  LOG     <- F
  yss     <- sd(as.vector(y1))/mean(y1)
  
  if(type == 'CA'){
    wpos <- length( which(y1 > 0) )
    nPerBin <-  max( c(10,wpos/15) )
  }
  if(type %in% c('PA', 'CAT')){
    breaks  <- c(-.05,.05,.95,1.05)
    wide    <- rep(.08,4)
    nPerBin <- NULL
    ylimit  <- c(0,1)
    xlimit <- c(-.1,1.1)
  } 
  if(type == 'OC'){
    breaks  <- seq(-.5,max(y1) + .5,by=1)
    wide    <- 1/max(y1)
    nPerBin <- NULL
    ylimit  <- range(yp,na.rm=T)
    xlimit  <- c( min(floor(y1)), max(ceiling(y1)) )
  } 
  if(type == 'DA')MEDIAN <- F
  if(type %in% c('DA','CA')){
    if(yss > 3){
      xlimit <- range(y1)
      xlimit[2] <- xlimit[2] + 1
      LOG <- T
    }
  }
  if(type %in% c('FC','CC')){
 #   xlimit <- range(y1)
    MEDIAN <- F
    nPerBin <- round( n*nk/50,0 )
 #   if(type  == 'CC')LOG <- T
  } 
  if(type == 'CC')xlimit[2] <- xlimit[2] + 1
  
  if( !is.null(censm) ){
    
    cc  <- censm$partition
    vlines  <- numeric(0)
    breaks  <- NULL
    nPerBin <- n*nk/15
    xlimit  <- range(y1)
    ylimit  <- quantile(yp,c(.01,.99),na.rm=T)
    
    if(ncol(cc) > 1){
      cm     <- unique( as.vector(cc[-1,]) )
      vlines <- cm[is.finite(cm)]
      breaks <- vlines
      nbin   <- nPerBin <- NULL
      uncens <- cbind(cc[3,-ncol(cc)],cc[2,-1])
      wu     <- which( uncens[,1] != uncens[,2] )
      for(m in wu){
        sm <- seq(uncens[m,1],uncens[m,2],length=round(10/length(wu),0))
        if(type == 'DA') sm <- c(uncens[m,1]:uncens[m,2])
        breaks <- c(breaks,sm)
      }
      if(max(cc[1,]) < Inf){
        breaks <- c(breaks, seq(max(breaks),(max(y1)+1),length=12) )
      } else {
        breaks <- c(breaks,max(y1) + 1)
      }
      breaks <- sort( unique(breaks) )
    }
  }
  
  if(LOG){
    xlimit[1] <- ylimit[1] <- 1
    w0     <- which(y1 == 0)
    y1[w0] <- ylimit[1]
    w0     <- which(yp == 0)
    yp[w0] <- ylimit[1]
    nPerBin <- nPerBin/4
    ylimit[1] <- xlimit[1] <- 1
    ylimit[2] <- max(yp)
  }
  
  list( y1 = y1, yp = yp, nbin=nbin, nPerBin=nPerBin, vlines=vlines,
        xlimit=xlimit,ylimit=ylimit,breaks=breaks,wide=wide,LOG=LOG,
        POINTS=F,MEDIAN=MEDIAN )
}
.gjamPredictTraits <- function(w,specByTrait,traitTypes){
  
  M  <- nrow(specByTrait)
  tn <- rownames(specByTrait)
  
  ww <- w
  ww[ww < 0] <- 0
  
  tt <- ww%*%t(specByTrait)
 # wf <- grep('FC',traitTypes)
 # if(length(wf) > 0){
 #   w0 <- which(tt[,wf] < 0)
 #   tt[tt[,wf] < 0,wf] <- 0
 #   tsum <- colSums(tt)
 #   tt   <- sweep(tt,1,tsum,'/')
 # }
  tt
}


.initW <- function(tw, x, yy, minw = -ncol(yy), cat=F){
  
  # initialize w for y = 0
  
  X <- x
  X[,-1] <- jitter(X[,-1],factor=1)
  
  XX  <- crossprod(X)
  IXX <- solve(XX)

    for(j in 1:50){
      
      bb <- IXX%*%crossprod(X,tw)
      muw <- X%*%bb
      
      tw[yy == 0] <- muw[yy == 0]    #neg values 
      tw[yy == 0 & tw > 0] <- 0      #no bigger than zero
    }
    tw[tw < minw] <- minw
 # }
  tw
}


.gjamSetup <- function(typeNames, x, y, breakList=NULL, holdoutN, holdoutIndex,
                       censor=NULL, effort=NULL, maxBreaks=100){
   
  Q <- ncol(x)
  n <- nrow(y)
  S <- ncol(y)
  
  effMat <- effort$values

  tmp <- .gjamGetTypes(typeNames)
  typeFull <- tmp$typeFull
  typeCols <- tmp$typeCols
  allTypes <- unique(typeCols)
  
  cuts <- cutLo <- cutHi <- numeric(0)
  minOrd <- maxOrd <- breakMat <- numeric(0)
  
  ordCols  <- which(typeNames == 'OC')
  disCols  <- which(typeNames == 'DA')
  compCols <- which(typeNames == 'CC')
  corCols  <- which(typeNames %in% c('PA','OC','CAT'))
  catCols  <- which(typeNames == c('CAT'))
  
  CCgroups <- attr(typeNames,'CCgroups')
  if(length(CCgroups) == 0)CCgroups <- rep(0,S)
  ngroup <- max(CCgroups)
  
  FCgroups <- attr(typeNames,'FCgroups')
  if(length(FCgroups) == 0)FCgroups <- rep(0,S)
  fgroup <- max(FCgroups)
  
  CATgroups <- attr(typeNames,'CATgroups')
  if(length(CATgroups) == 0)CATgroups <- rep(0,S)
  cgroup <- max(CATgroups)
  
  wo <- grep('others',colnames(y))
  if(length(wo) > 0){
    colnames(y)[wo] <- .replaceString(colnames(y)[wo],'others','other')
  }
  
  other <- grep('other',colnames(y))
  
  colnames(y) <- .replaceString(colnames(y),now=' ',new='')
  colnames(y) <- .replaceString(colnames(y),now='_',new='')
  colnames(x) <- .replaceString(colnames(x),now=' ',new='')
  
  w  <- y 
  if(!is.null(effort))w <- w/effort$values
  
  maxy <- apply(w,2,max,na.rm=T)
  maxy[maxy < 0] <- -maxy[maxy < 0]
  maxy <- matrix(maxy, n, S, byrow=T)
  
  z  <- w*0
  z[y == 0] <- 1
  z[y > 0]  <- 2
  plo <- phi <- y*0
  plo[z == 1] <- -2*maxy[z == 1]
  phi[z == 2] <- 2*maxy[z == 2]
  
  censorCA <- censorDA <- numeric(0) # CA and DA values to be sampled
  sampleW  <- y*0
  sampleW[is.na(sampleW)] <- 1
  
  for(k in allTypes){
    
    wk <- which(typeCols == k)
    nk <- length(wk)
    
    if( typeFull[wk[1]] == 'presenceAbsence' ){       
      
      sampleW[,wk] <- 1
      plo[,wk][z[,wk] == 1] <- -10
      phi[,wk][z[,wk] == 2] <- 10
      
      w[,wk] <- .tnorm(nk*n,plo[,wk],phi[,wk],0,1)
      br <- c(-Inf,0,Inf)
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('PA',wk,sep='-')
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat     <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'continuous' ){      
      
      sampleW[,wk] <- 0
      z[,wk]   <- 1
      
      phi[,wk] <- 5*maxy[,wk]
      plo[,wk] <- -phi[,wk]
      
      br <- c(-Inf,-Inf,Inf)
      br  <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('CON',wk,sep='-')
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'contAbun' ){       
      
      phi[,wk] <- 5*maxy[,wk]
      plo[,wk] <- -phi[,wk]
      
      w[,wk] <- .initW(w[,wk], x, y[,wk], minw = -max(y[,wk],na.rm=T)*5)
      
      br <- c(-Inf,0,Inf)
      br  <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('CA',wk,sep='-')
      
      sampleW[,wk][y[,wk] == 0] <- 1
      
      if( !is.null(censor) & 'CA' %in% names(censor) ){
        
        wc     <- which(names(censor) == 'CA')
        bc     <- censorCA <- numeric(0)
        
        for(m in wc){
          
          wm     <- censor[[m]]$columns
          cp     <- censor[[m]]$partition
          for(ii in 1:ncol(cp)){
            wmm <- which(y[,wm] == cp[1,ii] | (y[,wm] > cp[2,ii] & y[,wm] < cp[ii]) )
            plo[,wm][wmm] <- cp[2,ii]
            phi[,wm][wmm] <- cp[3,ii]
          }
          
          tmp    <- .gjamCensorSetup(y,w,z,plo,phi,wm,censorMat=
                                       censor[[m]]$partition)
     #     w[,wm] <- tmp$w[,wm]
          z[,wm] <- tmp$z[,wm]
     #     plo[,wm] <- tmp$plo[,wm]
     #     phi[,wm] <- tmp$phi[,wm]
          censorCA <- c(censorCA,tmp$censValue)
          bt       <- tmp$breakMat
          colnames(bt) <- as.character(c(1:ncol(bt)))
          rownames(bt) <- paste('CA',wm,sep='-')
          
          bc <- .appendMatrix(bc,bt,SORT=T,asNumbers=T)
        }
        
        mm <- match(rownames(bc),rownames(br))
        
        if(is.na(min(mm)))stop('error in censor list, check for conflicts')
        
        bb <- br[-mm,]
        tmp <- .appendMatrix(bc,bb,SORT=T,asNumbers=T)
        o   <- as.numeric( matrix( unlist(strsplit(rownames(tmp),'-')),
                                   ncol=2,byrow=T)[,2] ) 
        br <- tmp[order(o),]
      }
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'discAbun' ){
      
      plo[,wk] <- (y[,wk] - .5)/effMat[,wk]
      phi[,wk] <- (y[,wk] + .5)/effMat[,wk]
      plo[,wk][y[,wk] == 0] <- -5*maxy[,wk][y[,wk] == 0] 
      phi[,wk][y[,wk] == maxy[,wk]] <- 5*maxy[,wk][y[,wk] == maxy[,wk]]
      
      
      sampleW[,wk] <- 1
      
      disCols <- wk
      
      z[,wk] <- y[,wk] + 1
      w[,wk] <- .tnorm(nk*n,plo[,wk],phi[,wk],w[,wk],1)
      
      n <- nrow(y)
      S <- ncol(y)
      
      br <- c(-Inf,seq(0,(max(y[,wk])-1)),Inf)
      if(length(br) > maxBreaks){
        #     warning('breaks created')
        br <- c(br[1:maxBreaks],Inf)
      }
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('DA',wk,sep='-')
      
      if( !is.null(censor) & 'DA' %in% names(censor) ){
        
        wc     <- which(names(censor) == 'DA')
        bc     <- censorDA <- numeric(0)
        
        for(m in wc){
          wm     <- censor[[m]]$columns
          tmp    <- .gjamCensorSetup(y,w,z,plo,phi,wm,
                                     censorMat=censor[[m]]$partition)
          w[,wm] <- tmp$w[,wm]
          z[,wm] <- tmp$z[,wm]
          plo[,wm] <- tmp$plo[,wm]
          phi[,wm] <- tmp$phi[,wm]
          censorDA <- c(censorDA,tmp$censValue)
          bt       <- tmp$breakMat
          colnames(bt) <- as.character(c(1:ncol(bt)))
          rownames(bt) <- paste('DA',wm,sep='-')
          
          bc <- .appendMatrix(bc,bt,SORT=T,asNumbers=T)
        }
        mm <- match(rownames(bc),rownames(br))
        
        bb <- br[-mm,]
        tmp <- .appendMatrix(bc,bb,SORT=T,asNumbers=T)
        o   <- as.numeric( matrix( unlist(strsplit(rownames(tmp),'-')),
                                   ncol=2,byrow=T)[,2] ) 
        br <- tmp[order(o),]
      }
      
      
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'fracComp' ){
      
      wss <- which(y[,wk] == 0 | y[,wk] == 1)
      sampleW[,wk][wss] <- 1
      
      for(i in 1:fgroup){
        
        if(fgroup == 1){
          wki <- wk
        } else {
          wki <- which(typeCols == k & FCgroups == i)
        }
        nki <- length(wki)
        yki  <- y[,wki]
        
        lo <- plo[,wki]
        hi <- phi[,wki]
        lo[lo < -2/S] <- -20/S
        hi[hi > 3]  <- 3
        plo[,wki] <- lo
        phi[,wki] <- hi
        
        w[,wki] <- .initW(w[,wki], x, yki, minw = -20/S)
      }
      
      br <- c(-1,0,1)
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('FC',wk,sep='-')
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] %in% c('countComp','categorical')){
      
      sampleW[,wk] <- 1
      
      ntt <- ngroup
      if(typeFull[wk[1]] == 'categorical')ntt <- cgroup
      
      for(i in 1:ntt){
        
        if(ntt == 1){
          wki <- wk
        } else {
          wki <- which( typeCols == k )
          wki <- wki[ CCgroups[wki] == i | CATgroups[wki] == i ]
        }
        nki <- length(wki)
        yki  <- y[,wki]
        
        if( wki[1] %in% catCols ){
          lo  <- hi <- yki*0
          lo[yki == 0] <- -100
          hi[yki == 0] <- 0
          hi[yki == 1] <- 100
          mu <- yki*0
          mu[lo == 0] <- 20
          mu[hi == 0] <- -20
        } else {
          ee <- rowSums(yki)  + 1
          lo <- (yki - .5)/ee          
          hi <- (yki + .5)/ee
          lo[lo < 0]  <- -20/S
          mu <- (lo + hi)/2
        }
        
        z[,wki] <- yki + 1
        
        plo[,wki] <- lo
        phi[,wki] <- hi
        
        tmp <- matrix( .tnorm(nki*n,as.vector(lo),as.vector(hi), mu,sig=5),n,nki )
        
        tt <- tmp
        if(!wki[1] %in% catCols){
          tt[tt < 0] <- 0
          tsum <- rowSums(tt)
          tt   <- sweep(tt,1,tsum,'/')
          tt[tmp < 0] <- tmp[tmp < 0]
        }
        
   #     w[,wki] <- .initW(tt,x,y[,wki], minw = -100, cat=T)
        w[,wki] <- tt
      }
      
      br <- c(-1,0,1)
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('CC',wk,sep='-')
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'ordinal' ){
      
      nc <- apply(y[,wk,drop=F],2,max)
      
      sampleW[,wk] <- 1
      
      # more than one obs needed in last cell to estimate partition
      ii  <- list(spec = as.vector(matrix(c(1:nk),n,nk,byrow=T)), 
                  ss = as.vector(y[,wk,drop=F]))
      ctmp <- .byIndex(as.vector(y[,wk,drop=F])*0+1,ii,sum)
      
      ncc <- nc + 1
      if(max(ncc) > ncol(ctmp))ncc <- nc
      
      maxOne <- which(ctmp[ cbind(1:nk,ncc) ] == 1)
      
      if(length(maxOne) > 0){

        for(m in 1:length(maxOne)){
          mc <- wk[maxOne[m]]
          y[y[,mc] == nc[maxOne[m]],mc] <- nc[maxOne[m]] - 1
        }
        nc <- apply(y[,wk,drop=F],2,max)
      }
      
      ncut <- max(y[,wk,drop=F])
      crow <- c(0:ncut)
      cuts <- t( matrix(crow,(ncut+1),nk) )
      cuts[ cbind((1:nk),nc+1) ] <- Inf
      
      call <- t( apply(cuts,1,cumsum) )
      cuts[call == Inf] <- Inf
      cuts <- cbind(-Inf,cuts)
      if(!is.matrix(cuts))cuts <- matrix(cuts,1)
    
      tmp   <- .gjamGetCuts(y + 1,wk)
      cutLo <- tmp$cutLo
      cutHi <- tmp$cutHi
      
      ss   <- seq(0,(nk-1)*n,by=n)
      wh <- as.vector( outer(holdoutIndex,ss,'+') )
      c1 <- cutLo
      if(length(wh) > 0)c1 <- cutLo[-wh,]

      otab <- .byIndex(c1[,1]*0 + 1,INDICES=list('i'=c1[,1],
                                                 'j'=c1[,2]),sum,coerce=T)
      oo <- cbind(0,t( apply(otab,1,cumsum) ))
      wo <- which(oo == 0,arr.ind=T)
      wo[,2] <- as.numeric(colnames(otab))[wo[,2]]
      minOrd <- .byIndex(wo[,2],wo[,1],max)
      
      oo <- cbind(0,t( apply( t(apply(otab,1,rev)),1,cumsum) ))
      wo <- which(oo == 0,arr.ind=T)
      maxOrd <- ncut - .byIndex(wo[,2],wo[,1],max) + 2
      
      plo[,wk] <- cuts[cutLo]
      phi[,wk] <- cuts[cutHi]
      
      z[,wk] <- y[,wk] + 1
      w[,wk] <- matrix( .tnorm(nk*n,plo[,wk],phi[,wk],y[,wk],1),n,nk )
      
      colnames(cuts) <- c(1:ncol(cuts))
      rownames(cuts) <- paste('OC',wk,sep='-')
      rownames(cuts) <- paste(colnames(y)[wk],rownames(cuts),sep='_')
      breakMat <- .appendMatrix(breakMat,cuts,SORT=T,asNumbers=T)
    }
  }
  
  sord <- matrix( unlist(strsplit(rownames(breakMat),'_')),ncol=2,byrow=T)[,1]
  yord <- match(colnames(y),sord)
  breakMat <- breakMat[yord,]
  
  wCols <- which(colSums(sampleW) > 0)
  wRows <- which(rowSums(sampleW) > 0)
  
  attr(sampleW,'type') <- 'cols'
  attr(sampleW,'index') <- wCols
  if( sum(sampleW) == 0)attr(sampleW,'type') <- 'none'
  if( sum(sampleW) > 0 & (length(wRows) < length(wCols)) ){
    attr(sampleW,'type') <- 'rows'
    attr(sampleW,'index') <- wRows
  }
  
  ii <- list(spec = as.vector(matrix(c(1:S),n,S,byrow=T)), 
             discrete_class = as.vector(z))
  classBySpec <- .byIndex(as.vector(z)*0+1,ii,sum)
  rownames(classBySpec) <- colnames(y)
  
  ncc <- min(20,ncol(classBySpec))
  nrr <- min(20,nrow(classBySpec))
 # print( classBySpec[1:nrr,1:ncc] )
  
    list(w = w, z = z, y = y, other = other, cuts = cuts, 
         cutLo = cutLo, cutHi = cutHi, 
         plo = plo, phi = phi, ordCols=ordCols, disCols = disCols, 
         compCols = compCols, corCols = corCols,
         classBySpec = classBySpec, breakMat = breakMat, 
         minOrd = minOrd, maxOrd = maxOrd, sampleW = sampleW,
         censorCA = censorCA, censorDA = censorDA )
 }

.gjamTheta2cuts <- function(tg,ss){
  nc   <- ncol(tg)
  sr    <- nrow(ss)
  if(length(ss) == 1)return( tg/sqrt(ss) )
  tg/matrix( sqrt(diag(ss)),sr,nc)
}

.gjamTrueVest <- function(chains,true,typeCode,allTypes,xlim=NULL,ylim=NULL,
                          label=NULL,colors=NULL,add=F,legend=T){
  
  true   <- as.vector(true)
  ntypes <- length(allTypes)
  
  if(is.null(ylim))ylim <- range(chains,na.rm=T)
  if(is.null(xlim))xlim <- range(true,na.rm=T)
  
  if(!is.matrix(chains)){
    chains <- matrix(chains,ncol=1)
    bCoeffTable <- c(mean(chains),sd(chains),quantile(chains,c(.025,.975)),true)
    bCoeffTable <- matrix(bCoeffTable,1)
  } else {
    bCoeffTable <- .processPars(chains,xtrue=true )
  }
  
  if(is.null(colors)){
    colors <- 1
    if(ntypes > 1)colors <- typeCode
  }
  if(length(colors) == 1) colors <- rep(colors,ntypes)
  
  .predVsObs(true,p=chains,xlab='true',xlim=xlim,ylim=ylim,ylab='estimated',
            colors=colors,add=add)
  
  if(ntypes > 1 & legend)legend('topleft',allTypes,text.col=colors,bty='n')
  if(!is.null(label)).plotLabel(label,above=T)
  
  invisible( bCoeffTable )
}
.gjamUpdateBetaNoPrior <- function(WIX,IXX,sg,...){
  
  bg <- matrix( .rMVN(1,as.vector(WIX),kronecker(sg,IXX)),nrow(IXX),ncol(WIX) )
  
  list(bg = bg)
}

.conditionalMVNRcpp <- function(xx, mu, sigma, cdex, p=ncol(mu)){  
  # xx, mu are matrices
  
  # cdex conditional for these variables
  # gdex condition on these variables
  
  if(ncol(xx) != ncol(sigma))stop('ncol(xx) != ncol(sigma)')
  if(ncol(mu) != ncol(sigma))stop('ncol(mu) != ncol(sigma)')
  if(max(cdex) > ncol(mu))stop('max(cdex) > ncol(mu)')
  
  # if(!length(xx) > nrow(sigma))return( conditionalMVN(xx,mu,sigma,cdex) )
  
  gdex <- (1:p)[-cdex] - 1
  cdex <- cdex - 1
  .conditionalMVNRcppCpp(cdex, gdex, xx, mu, sigma) 
}

.byRcpp <- function(x, i, j, summat=matrix(0,max(i),max(j)), 
                   totmat=summat, fun='mean'){  #
  
  nn <- length(x)
  if( nn != length(i) | nn != length(j) )stop('vectors unequal in byFunctionRcpp')
  if( nrow(summat) < max(i) | ncol(summat) < max(j) )stop('matrix too small')
  
  ww <- which(is.na(x))
  if(length(ww) > 0){
    x <- x[-ww]
    i <- i[-ww]
    j <- j[-ww]
  }
  
  frommat <- cbind(i,j,x)
  
  nr  <- nrow(frommat)
  maxmat <- summat*0 - Inf
  minmat <- summat*0 + Inf
  
  tmp <- .byRccpCpp(nr, frommat, totmat, summat, minmat, maxmat)
  
  if(fun == 'sum')return(tmp$sum)
  if(fun == 'mean'){
    mu <- tmp$sum/tmp$total
    mu[is.na(mu)] <- 0
    return(mu)
  }
  if(fun == 'min'){
    return( tmp$min )
  }
  tmp$max
}
.tnormMVNmatrixRcpp <- function(avec, muvec, smat, 
                                lo=matrix(-1000,nrow(muvec),ncol(muvec)), 
                                hi=matrix(1000,nrow(muvec),ncol(muvec)),
                                whichSample = c(1:nrow(smat))){
  
  #lo, hi must be same dimensions as muvec,avec
  
  lo[lo < -1000] <- -1000
  hi[hi > 1000]  <- 1000
  
  if(max(whichSample) > length(muvec))stop('whichSample outside length(muvec)')
  
  r <- avec
  a <- .tnormMVNmatrixRcppCpp(avec, muvec, smat, lo, hi, whichSample, 
                              idxALL = c(0:(nrow(smat)-1)) )  
  r[,whichSample] <- a[,whichSample]
  r
}

.whichFactor <- function(dframe){
  
  if(!is.data.frame(dframe))return(character(0))
  
  tmp <- model.frame(data = dframe)
  ym <- attr( attributes(tmp)$terms, 'dataClasses' )
  
  which(ym == 'factor')
}


.xpredSetup <- function(w, x, bg, isNonLinX, isFactor, intMat, standMat, standMu,
                        factorList, notOther, notStandard){
  
  # initialize predicted X
  
  xpred  <- x
  n      <- nrow(x)
  Q      <- ncol(x)
  xnames <- colnames(x)
  SO     <- length(notOther)
  
  px <- 1:Q
  if(length(isNonLinX) > 0)px <- px[-isNonLinX]
  px <- px[!xnames[px] %in% isFactor]
  px <- px[px != 1]
  
  ii <- grep(':',xnames,fixed=T)
  i2 <- grep('^2',xnames,fixed=T)
  
  qx <- c( 1, ii, i2)
  qx <- c(1:Q)[-qx]
  bx <- bg[drop=F,qx,notOther]
  cx <- crossprod(t(bx))
  if(length(cx) == 1){
    cx <- 1/(cx*1.001)
  } else {
    cx <- cx + diag(diag(cx)*.001)
    cx <- solve(cx)
  }
  
  xx <- (w[,notOther] - matrix(bg[1,notOther],n,SO,byrow=T))%*%t(bx)%*%cx
  colnames(xx) <- xnames[qx]
  scol      <- colnames(xx)[!colnames(xx) %in% notStandard]
  xx[,scol] <- sweep(xx[,scol,drop=F],2,colMeans(xx[,scol,drop=F]),'-')
  xx[,scol] <- sweep(xx[,scol,drop=F],2,apply(xx[,scol,drop=F],2,sd),'/')
  xpred[,qx] <- xx
  
  # x is centered/standardized
#  xx <- standMu[,1] + 2*standMat[,1]
#  xx <- matrix(xx,n,Q,byrow=T)
#  xpred[xpred > xx] <- xx[xpred > xx]
  
#  xx <- standMu[,1] - 2*standMat[,1]
#  xx <- matrix(xx,n,Q,byrow=T)
#  xpred[xpred < xx] <- xx[xpred < xx]
  
  xpred[xpred < -4] <- -4
  xpred[xpred > 4] <- 4
  
  if(length(intMat) > 0){
    for(k in 1:nrow(intMat)){
      xpred[,intMat[k,1]] <- xpred[,intMat[k,2]]*xpred[,intMat[k,3]]
    }
  }
  
  linFactor <- numeric(0)
  
  if(length(isFactor) > 0){
    xpred[,isFactor] <- 0
    
    for(k in 1:length(factorList)){
      kf  <- lf <- factorList[[k]]
      
      if( !is.null(isNonLinX) ){
        xin <- xnames[isNonLinX]
        lf  <- kf[!kf %in% xin]
      }
      if(length(lf) == 0)next
      lf  <- match(lf,xnames)
      ww  <- which(is.finite(lf))
      
      wt  <- colSums(x[,c(1,lf)])   #random, but weighted by prevalence
      wt  <- wt/sum(wt)
      sk  <- sample(c(1,lf),n, replace=T, prob=wt)
      xpred[ cbind(c(1:n),sk) ] <- 1
      
      if(length(ww) == 0)next
      lf <- c(1,lf)   # intercept is reference
      linFactor <- append(linFactor, list(lf))
    }
  }
  
  list(linFactor = linFactor, xpred = xpred, px = px)
}

.blockDiag <- function(mat1,mat2){
  
  #creates block diagional
  
  if(length(mat1) == 0)return(mat2)
  
  namesc <- c(colnames(mat1),colnames(mat2))
  namesr <- c(rownames(mat1),rownames(mat2))
  
  nr1 <- nrow(mat1)
  nr2 <- nrow(mat2)
  nc1 <- ncol(mat1)
  nc2 <- ncol(mat2)
  nr  <- nr1 + nr2
  nc  <- nc1 + nc2
  
  new <- matrix(0,nr,nc)
  new[ 1:nr1, 1:nc1 ] <- mat1
  new[ (nr1+1):nr, (nc1+1):nc ] <- mat2
  colnames(new) <- namesc
  rownames(new) <- namesr
  new
}

.getContrasts <- function(facK, fnames){
  
  # D - x to z
  # L - beta to alpha
  # facK - name of factor
  # fnames - character of factor levels
  
  ff <- paste(facK,fnames,sep='')
  
  Q  <- length(fnames)
  cc <- diag(Q)
  cc[1,] <- -1
  dd <- cc
  dd[1] <- 1
  cc[,1] <- 1
  colnames(cc) <- colnames(dd) <- c('intercept',ff[-1])
  rownames(cc) <- rownames(dd) <- fnames
  L <- t( solve(cc) )
  list(C = cc, D = dd, L = L)
}

.getUnstandX <- function(xx, standRows, xmu, xsd, intMat){
  # coefficients to unstandard scale
  
  n   <- nrow(xx)
  nsr <- length(standRows)
  
  xsm <- matrix(xsd,n,nsr,byrow=T)
  
  unstand  <- matrix(xmu,n,nsr,byrow=T) + xx[,standRows,drop=F]*xsm
  xUnstand <- xx
  xUnstand[,standRows] <- unstand
  if(length(intMat) > 0){
    for(j in 1:nrow(intMat)){
      xUnstand[,intMat[j,1]] <- xx[,intMat[j,2]] * xx[,intMat[j,3]] 
    }
  }
 
  S2U <- ginv(crossprod(xUnstand))%*%t(xUnstand)
  list(xu = xUnstand, S2U = S2U)
}

.getHoldLoHi <- function(yh, wh, pl, ph, eff, ymax, typeNames, cutg, ordCols){
  
  # update plo, phi for holdouts, yh is prediction
  
  allTypes <- unique(typeNames)
  
  for(k in 1:length(allTypes)){
    
    tk <- allTypes[k]
    wk <- which(typeNames == tk)
    
    if(tk == 'CON')next
    
    #CAT
    
  #  if(tk == 'OC'){
  #    tmp   <- .gjamGetCuts(yh+1,ordCols)  #holdout pred y cuts
  #    pl[,ordCols] <- cutg[tmp$cutLo]
  #    ph[,ordCols] <- cutg[tmp$cutHi]
  #  }
    if(tk == 'PA'){
      pl[,wk][yh[,wk] == 0] <- -10
      pl[,wk][yh[,wk] == 1] <- 0
      ph[,wk][yh[,wk] == 0] <- 0
      ph[,wk][yh[,wk] == 1] <- 10
    }
    if(tk == 'CA'){
      ym <- max(ymax[wk])
      pl[,wk][yh[,wk] == 0] <- -5*ym
      pl[,wk][yh[,wk] > 0]  <- 0
      ph[,wk][yh[,wk] == 0] <- 0
      ph[,wk][yh[,wk] > 0] <- 5*ym
    }
    if(tk == 'DA'){
      ym <- max(ymax[wk])
      ee <- 1
      if(!is.null(eff))ee <- eff[,wk]
      pl[,wk] <- (yh[,wk] - .5)/ee
      ph[,wk] <- (yh[,wk] + .5)/ee
      pl[,wk][yh[,wk] == 0] <- -5*ym
      pl[,wk][yh[,wk] == ym] <- 5*ym
    }
    if(tk == 'FC'){
      pl[,wk][yh[,wk] == 0] <- -5
      pl[,wk][yh[,wk] > 0]  <- 0
      pl[,wk][yh[,wk] > 1]  <- 1
      ph[,wk][yh[,wk] == 0] <- 0
      ph[,wk][yh[,wk] > 0]  <- 1
      ph[,wk][yh[,wk] == 1] <- 5
    }
    if(tk == 'CC'){
      ym <- rowSums(yh[,wk,drop=F])
      ee <- matrix(ym,nrow(yh),length(wk))
      pl[,wk] <- (yh[,wk] - .5)/ee
      ph[,wk] <- (yh[,wk] + .5)/ee
      pl[,wk][yh[,wk] == 0] <- -5
      pl[,wk][yh[,wk] == ee] <- 5
    }
  }
  
  list(pl = pl, ph = ph)
}

.gjam <- function(formula, xdata, ydata, modelList) UseMethod(".gjam")

.gjam.default <- function(formula, xdata, ydata, modelList){
  
  holdoutN      <-  0
  holdoutIndex  <- numeric(0)
  modelSummary  <- betaPrior  <- traitList <- effort <- NULL
  specByTrait   <- traitTypes <- breakList <- notStandard <- NULL
  censor <- censorCA <- censorDA <- CCgroups <- FCgroups <- NULL
  y0 <- N  <- r <- otherpar <- pg <- NULL
  ng     <- 2000
  burnin <- 500
  ZEROINFL <- REDUCT <- TRAITS <- F
  PREDICTX <- T
  testCC   <- F     #test w's in CC model
  
  ematAlpha <- .5

  alpha.DP <- ncol(ydata)          # large values give more variation
  
  if(alpha.DP == 1)stop('multivariate model: at least 2 columns needed in ydata')
  
  for(k in 1:length(modelList))assign( names(modelList)[k], modelList[[k]] )
  if(!is.null(traitList)){
    TRAITS <- T
    for(k in 1:length(traitList))assign( names(traitList)[k], traitList[[k]] )
    
    stt <- .replaceString(colnames(specByTrait),'_','')
    colnames(specByTrait) <- stt
    colnames(plotByTrait) <- stt
    colnames(traitList$specByTrait) <- stt
    colnames(traitList$plotByTrait) <- stt
    modelList$traitList <- traitList
  }
  
  if(burnin >= ng)           stop( 'burnin must be > no. MCMC steps, ng' )
  if('censor' %in% names(modelList)){
    for(k in 1:length(censor)){
      if(!names(censor)[[k]] %in% c('CA','DA'))
        stop('censor name must be CA or DA')
      if( nrow(censor[[k]]$partition) != 3 )
        stop('censor matrix: 3 rows for value, lo, hi')
      rownames(censor[[k]]$partition) <- c('value','lo','hi')
    }
  }
  
  if(missing(xdata)) xdata <- environment(formula)
  
  
  S <- ncol(ydata)
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  if(length(typeNames) != S) 
    stop('typeNames must be one value or no. columns in y')
 # if( !all( !sapply(as.data.frame(ydata),is.factor) )  )
 #   stop('ydata cannot contain factors')
  if( !all( !sapply(as.data.frame(ydata),is.character) )  )
    stop('ydata cannot contain characters')
  
  if(TRAITS){
    if(!all( typeNames %in% c('CC','FC') ) )
      stop('trait prediction requires composition data (CC or FC)')
    if(nrow(plotByTrait) != nrow(ydata))
      stop('nrow(plotByTrait) must equal nrow(ydata)')
    if(ncol(plotByTrait) != length(traitTypes))
      stop('ncol(plotByTrait) must equal length(traitTypes)')
    if(ncol(plotByTrait) != length(traitTypes))
      stop('ncol(plotByTrait) must equal length(traitTypes)')
    ii <- identical(rownames(specByTrait),colnames(ydata))
    if(!ii){
      ww <- match(colnames(ydata),rownames(specByTrait) )
      if( is.finite(min(ww)) ){
        specByTrait <- specByTrait[ww,]
      } else {
        stop( 'rownames(specByTrait) must match colnames(ydata)' )
      }
    }
    if(typeNames[1] == 'CC')ydata <- round(ydata,0)
  }
  
  tmp <- .buildYdata(ydata, typeNames)
  y   <- tmp$y
  ydataNames <- tmp$ydataNames
  typeNames  <- tmp$typeNames
  CCgroups   <- tmp$CCgroups
  FCgroups   <- tmp$FCgroups
  CATgroups  <- tmp$CATgroups
  if(TRAITS) rownames(specByTrait) <- colnames(y)
    
  S <- ncol(y)
  n <- nrow(y)
  
  effMat <- y*0 + 1
  effMat[is.na(effMat)] <- 0
  if( is.null(effort)){
    effort <- list(columns = 1:S, values = effMat)
  } else {
    effMat[,effort$columns] <- effort$values
    effort$values <- effMat
  }
  effort$columns <- 1:S
  
  tmp      <- .gjamGetTypes(typeNames)
  typeCols <- tmp$typeCols
  typeFull <- tmp$typeFull
  typeCode <- tmp$TYPES[typeCols]
  allTypes <- sort(unique(typeCols))
  
  tmp <- .gjamXY(formula, xdata, y, typeNames, notStandard)
  x   <- tmp$x; y <- tmp$y; snames <- tmp$snames
  xdata  <- tmp$xdata
  xnames <- tmp$xnames
  isInt         <- tmp$isInt;        intMat <- tmp$intMat
  factorList    <- tmp$factorList; isFactor <- tmp$isFactor
  isNonLinX     <- tmp$isNonLinX
  designTable   <- tmp$designTable;  xscale <- tmp$xscale
  predXcols     <- tmp$predXcols
  standMat      <- tmp$standMat;    standMu <- tmp$standMu
  standRows     <- tmp$standRows;  facNames <- tmp$facNames
  contrast      <- tmp$contrast; xdataNames <- tmp$xdataNames
  
  if(length(factorList) == 0)factorList <- NULL

  modelSummary  <- append(modelSummary,
                          list(isFactor = isFactor, 
                               factorList = factorList, contrasts = contrast))
  
  modelList     <- append(modelList, list('formula' = tmp$formula,
                                          'notStandard' = notStandard))
  
  Q <- ncol(x)
  
  tmp <- .gjamMissingValues(x,y)
  xmiss  <- tmp$xmiss; xbound <- tmp$xbound; 
  ymiss  <- tmp$ymiss; missY <- tmp$missY
  xprior <- tmp$xprior
  yprior <- tmp$yprior
  nmiss  <- length(xmiss)
  mmiss  <- nrow(ymiss)
  if(nmiss > 0)x[xmiss] <- xprior
  if(mmiss > 0)y[ymiss] <- yprior
  
  npar  <- S*(S + 1)/2 
  nobs  <- S*n
  ratio <- 1/5
  Smax  <- floor( 2*n*ratio - 1 )
  Nmax  <- min( floor( c(S*n*ratio/3 , n/2)) )    # r  = 3
  
  OVERRIDE <- F
  if( 'REDUCT' %in% names(modelSummary) ){
    if( !modelSummary$REDUCT ) OVERRIDE <- T
  }
  if( !'reductList' %in% names(modelList) & S > min(Smax,100) & !OVERRIDE ){
    r <- max(c(3,Nmax/2))
    r <- floor(r)
    if(r > 8)r <- 8
    reductList <- list(r = r, N = Nmax, alpha.DP = alpha.DP)
    for(k in 1:length(reductList))assign( names(reductList)[k], reductList[[k]] )
    REDUCT <- T
  }
  if('reductList' %in% names(modelList) & !OVERRIDE){
    for(k in 1:length(reductList))assign( names(reductList)[k], reductList[[k]] )
    REDUCT <- T
  }
  
  
  tmp <- .gjamHoldoutSetup(holdoutIndex, holdoutN, n)
  holdoutIndex <- tmp$holdoutIndex; holdoutN <- tmp$holdoutN
  inSamples    <- tmp$inSamples; nIn <- tmp$nIn
  
  
  tmp <- .gjamSetup(typeNames, x, y, breakList, holdoutN, holdoutIndex,
                    censor=censor, effort=effort) 
  w <- tmp$w; z <- tmp$z; y <- tmp$y; other <- tmp$other; cuts <- tmp$cuts
  cutLo         <- tmp$cutLo; cutHi     <- tmp$cutHi
  plo <- tmp$plo; phi <- tmp$phi
  ordCols     <- tmp$ordCols; disCols <- tmp$disCols
  compCols    <- tmp$compCols 
  classBySpec <- tmp$classBySpec; breakMat <- tmp$breakMat
  minOrd      <- tmp$minOrd; maxOrd <- tmp$maxOrd; censorCA <- tmp$censorCA
  censorDA    <- tmp$censorDA; ncut <- ncol(cuts);  corCols <- tmp$corCols
  catCols     <- which(attr(typeNames,'CATgroups') > 0)
  sampleW     <- tmp$sampleW
  
  sampleW[censorCA] <- 1
  sampleW[censorDA] <- 1
  sampleWhold <- tgHold <- NULL
  wHold <- NULL
  wmax  <- apply(y/effMat,2,max,na.rm=T)
  pmin  <- -2*wmax
  
  wlo <- which(plo < pmin,arr.ind=T)
  if(length(wlo) > 0) plo[wlo] <- pmin[wlo[,2]]
  
  if(mmiss > 0){
    phi[ ymiss ] <- wmax[ ymiss[,2] ]
    plo[ ymiss ] <- pmin[ ymiss[,2] ]
    sampleW[ ymiss ] <- 1
  }
  
  ploHold <- phiHold <- NULL
  
  if(holdoutN > 0){
    sampleWhold <- sampleW[holdoutIndex,]  #to predict X
    sampleW[holdoutIndex,] <- 1
    tgHold  <- cuts
    wHold   <- w[drop=F,holdoutIndex,]
    ploHold <- plo[drop=F,holdoutIndex,]
    phiHold <- phi[drop=F,holdoutIndex,]
  }

  byCol <- byRow <- F
  if(attr(sampleW,'type') == 'cols')byCol <- T
  if(attr(sampleW,'type') == 'rows')byRow <- T
  indexW <- attr(sampleW,'index')
  
  notCorCols <- c(1:S)
  if(length(corCols) > 0)notCorCols <- notCorCols[-corCols]
  
  
  ############ beta
  updateBeta <- .gjamUpdateBetaNoPrior
  loBeta <- hiBeta <- NULL
  BPRIOR <- F
  
  if( !is.null(betaPrior) ){
    loBeta <- matrix(-10000,Q,S)
    hiBeta <- matrix(10000,Q,S)
    rownames(loBeta) <- rownames(hiBeta) <- xnames
    colnames(loBeta) <- colnames(hiBeta) <- snames
    
    wm <- match(rownames(betaPrior$lo),xnames)
    
    loBeta[wm,]     <- betaPrior$lo
    hiBeta[wm,]     <- betaPrior$hi
    updateBeta <- .gjamUpdateBetaPrior
    BPRIOR <- T
  }                 

  ############ 'other' columns
  sigmaDf  <- nIn - Q + S - 1
  sg <- diag(.1,S)
  SO <- S
  
  notOther <- c(1:S)
  sgOther  <- NULL
  if(length(other) > 0){                     
    notOther   <- notOther[!notOther %in% other]
    SO         <- length(notOther)
    sg[other,] <- sg[,other] <- 0
    sgOther    <- matrix( cbind(other,other),ncol=2 )
    sg[sgOther] <- .0001
  }
  modelSummary <- append(modelSummary,list(classBySpec = classBySpec))
  
  if(byCol){
    inw <- intersect( colnames(y)[indexW], colnames(y)[notOther] )
    indexW <- match(inw,colnames(y)[notOther])
  }
    
  loB <- hiB <- NULL
  if(BPRIOR){
    loB <- loBeta[,notOther]
    hiB <- hiBeta[,notOther]
  }
  
  ############ dimension reduction
    
  .param.fn <- .paramWrapper(REDUCT, inSamples, SS=length(notOther),
                             loB = loB, hiB = hiB, updateBeta)
  sigmaerror <- .1
  otherpar   <- list(S = S, Q = Q, sigmaerror = sigmaerror, Z = NA, K =rep(1,S),
                     sigmaDf = sigmaDf)
  ogibbs <- chiSum <- kgibbs <- sigErrGibbs <- rndEff <- NULL
  
  bg  <- alpha <- matrix(0,Q,S)
  
  
  yp <- y
  wmax <- ymax <- apply(y,2,max)
  wmax <- wmax/effMat

  if(REDUCT){
    message( paste('Dimension reduced from',S,'X',S,'to',N,'X',r,'responses') )
    SelPars           <- list()
    SelPars$p0.cols   <- 1                 # always have intercept
    SelPars$GammaList <- list()
    nbetasel          <- (Q - SelPars$p0.cols)*SO
    SelPars$GammaList$chi <- rbinom(nbetasel,1,0.5)
    SelPars$GammaList$tau <- rep(1,nbetasel)
    SelPars$GammaList$nu0 <- 0.0001
    SelPars$GammaList$omega    <- rep(0.5,nbetasel)
    SelPars$GammaList$pars.tau <- c(0.5,0.5)
    otherpar$SelPars <- SelPars
    
    otherpar$N <- N
    otherpar$r <- r
    otherpar$Z <- .rmvnormArma(N,rep(0,r),1/S*diag(r))
    otherpar$D <- .riwish(df = (2 + r + N), 
                          S = (crossprod(otherpar$Z) +
                                 2*2*diag(rgamma(r,shape=1,rate=0.001))))
    otherpar$K <- sample(1:N,length(notOther),replace=T)
    otherpar$sigmaerror <- 0.1
    otherpar$alpha.DP <- alpha.DP
    otherpar$pvec     <- .sample.p(N=N, avec=rep(alpha.DP/N,(N-1)),
                                   bvec=((N-1):1)*alpha.DP/N, K=otherpar$K)
    kgibbs <- matrix(1,ng,S)
    sgibbs <- matrix(0,ng, N*r)
    nnames <- paste('N',1:N,sep='-')
    rnames <- paste('r',1:r,sep='-')
    colnames(sgibbs) <- .multivarChainNames(nnames,rnames)
    sigErrGibbs <- rep(0,ng)      

    bi       <- bg*0
    bi[1,]   <- 1
    bi[,other] <- 1
    selIndex <- which(bi == 0)
  } else {
    Kindex <- which(as.vector(lower.tri(diag(S),diag=T)))
    nK     <- length(Kindex)
    sgibbs <- matrix(0,ng,nK)
    colnames(sgibbs) <- .multivarChainNames(snames,snames)[Kindex] # half matrix
  }

  ############ parameters

  out <- .param.fn(x, beta = bg[,notOther], Y = w[,notOther], otherpar)
  
  sg[other,] <- sg[,other] <- 0
  sg[notOther,notOther]    <- out$sg
  diag(sg)[other]          <- 1
 
  bg[,notOther] <- alpha[,notOther] <- out$bg 
  rownames(bg)  <- xnames
  otherpar      <- out$otherpar
  rownames(sg)  <- colnames(sg) <- colnames(bg) <- snames
  colnames(x)   <- xnames
  
  ############  traits
  if(length(specByTrait) > 0){
    specTrait <- specByTrait[colnames(y),]
    tnames    <- colnames(specTrait)
    M         <- ncol(specTrait)
    specTrait <- t(specTrait)
    
    missTrait <- which(is.na(specTrait),arr.ind=T)
    if(length(missTrait) > 0){
      traitMeans <- rowMeans(specTrait,na.rm=T)
      specTrait[missTrait] <- traitMeans[missTrait[,2]]
      warning( paste('no. missing trait values:',nrow(missTrait)) )
    }
      
    
    agibbs <- matrix(0,ng,M*Q)
    mgibbs <- matrix(0,ng,M*M)
    tpred  <- tpred2 <- matrix(0,n,M)
    colnames(agibbs) <- .multivarChainNames(xnames,tnames)
    colnames(mgibbs) <- .multivarChainNames(tnames,tnames)
  }
  
  ############ ordinal data
  if('OC' %in% typeCode){
    tg       <- cutg <- cuts
    cnames   <- paste('C',1:ncut,sep='-')
    nor      <- length(ordCols)
    cgibbs   <- matrix(0,ng,(ncut-3)*nor)
    colnames(cgibbs) <- as.vector( outer(snames[ordCols],
                                         cnames[-c(1,2,ncut)],paste,sep='_') )
    tmp   <- .gjamGetCuts(z,ordCols)
    cutLo <- tmp$cutLo
    cutHi <- tmp$cutHi
    plo[,ordCols] <- tg[cutLo]                                        
    phi[,ordCols] <- tg[cutHi]
    lastOrd <- ncol(tg)
  }
  
  ############ setup w
  tmp <- .gjamGetTypes(typeNames)
  typeFull <- tmp$typeFull
  typeCols <- tmp$typeCols
  allTypes <- unique(typeCols)
  
  
  .updateW <- .wWrapper(REDUCT, n, S, effMat, corCols, typeNames, 
                        typeFull, typeCols, 
                        allTypes, holdoutN, holdoutIndex, censor, 
                        censorCA, censorDA, notOther, sampleW, byRow, byCol,
                        indexW, ploHold, phiHold, sampleWhold)
  ycount <- rowSums(y)
  if('CC' %in% typeCode)ycount <- rowSums(y[,compCols])
  
  rndEff <- w*0

  ############ X prediction
  tmp <- .xpredSetup(w, x, bg, isNonLinX, isFactor, intMat, standMat, standMu,
                                factorList, notOther, notStandard) 
  linFactor <- tmp$linFactor
  xpred     <- tmp$xpred
  px        <- tmp$px
  nfact     <- length(factorList)
  priorXIV  <- diag(1e-5,ncol(x))
  priorX    <- colMeans(x)
  lox       <- apply(x,2 ,min)
  hix       <- apply(x,2,max)
  
  lox[isFactor] <- -2.5
  hix[isFactor] <- 2.5
  if(length(intMat) > 0){
    lox[intMat[,1]] <- -2.5
    hix[intMat[,1]] <- 2.5
  }

  ws        <- which(notStandard %in% xnames)
  if(length(ws) == 0){
    notStandard <- NULL
  } else {
    notStandard <- notStandard[ws]
    lox[notStandard] <- standMu[notStandard,1] - 3*standMat[notStandard,1]
    hix[notStandard] <- standMu[notStandard,1] + 3*standMat[notStandard,1]
  }
  
  ############  contrasts, predict F matrix
  
  q1 <- Q - 1
  fnames <- xnames
  findex <- character(0)
  
  if(nfact > 0){              # exclude main effects of factors
    findex <- sort( unique( unlist(factorList) ) )
    fnames <- fnames[!fnames %in% findex]
  }
    
  tmp <- diag(length(fnames))
  rownames(tmp) <- colnames(tmp) <- fnames
  if(length(tmp) < 2){
    eCont <- frow <- intercept <- numeric(0)
  } else {
    eCont <- tmp[drop=F,-1,]
    frow  <- rep(0,nrow(eCont))
    intercept <- rep(0,nrow(eCont))
  }
  dCont <- lCont <- eCont

  if(nfact > 0){
    
    for(k in 1:nfact){
      
      cm <- contrast[[k]]
      colnames(cm) <- factorList[[k]]
      
      facK <- names(factorList)[[k]]
      
      wx <- match(facK,colnames(xdata))
      
      fnames <- as.character( levels(xdata[[wx]]) ) 
      mm     <- .getContrasts(facK, fnames)
      D  <- mm$D                      # for Z <- x%*%D; 
      L  <- mm$L                      # for A <- L%*%bg; 
      C  <- mm$C                      # L <- solve(t(C)); C = solve(t(L))
      
      if(length(eCont) > 1){
        eCont <- .blockDiag(eCont,cm)
        dCont <- .blockDiag(dCont,D[,-1,drop=F])
        lCont <- .blockDiag(lCont,L[,-1,drop=F])
        ec    <- nrow(lCont)
        bc    <- ec - nrow(L) + 1
        lCont[bc:ec,1] <- L[,1]
        dCont[bc,1] <- -1
      } else {
        eCont <- cbind(0,cm)
        colnames(eCont)[1] <- 'intercept'
        dCont <- D
        lCont <- L
  #      colnames(dCont)[1] <- colnames(lCont)[1] <- 'intercept'
      }
      nr2   <- nrow(cm)
      nc2   <- ncol(cm)
      intercept <- c(intercept,rep(1,nr2))
      
      frow <- c(frow,rep(k,nr2))
    }
    
    eCont[,1] <- intercept
    q1 <- nrow(eCont)
    
    for(k in 1:nfact){
      
      tmp <- .replaceString(rownames(eCont),
                                        paste(names(factorList)[[k]],sep=''),'')
      wk <- which(nchar(tmp) > 0)
      rownames(eCont)[wk] <- tmp[wk]
    }
  }
  
  eCont <- eCont[drop=F,,xnames]
  
  dCont <- t(dCont[drop=F,,xnames])
  dCont[1,] <- abs(dCont[1,])
  lCont <- lCont[drop=F,,xnames]
  
  q1 <- nrow(eCont)
  fnames   <- rownames(eCont)
  facList2 <- factorList
  if(nfact > 0){
    for(j in 1:nfact){
      wj <- which(names(xdata) == names(factorList)[j])
      facList2[[j]] <- levels(xdata[[wj]])
    }
  }

          
  fmat <- matrix(0,q1,q1)
  colnames(fmat) <- rownames(fmat) <- fnames
  modelSummary <- append(modelSummary, list(dCont = dCont, eCont = eCont, lCont = lCont))
  
  findex <- match(findex,xnames)
  
  ############ E matrix
  emat <- matrix(0,S,S)
  colnames(emat) <- rownames(emat) <- snames
  lo <- hi <- lm <- hm <- ess <- emat
  
  ############ sp richness
  richness <- richFull <- NULL
  RICHNESS <- F
  
  notRichness <- which(!typeNames %in% c('CON','CAT','OC'))
  if(length(notRichness) > 0)RICHNESS  <- T
  
  wrich <- y*0 
  wrich[,notRichness] <- 1
  wrich[ymiss] <- 0
 
  covx <- cov(x)
  
  ############ sums
  predx  <- predx2 <- x*0
  yerror <- ypred  <- ypred2 <- wpred  <- wpred2 <- ymissPred <- ymissPred2 <- y*0
  sumDev <- 0   #for DIC
  sMean  <- sg*0
  ntot   <- 0
  
  ############ gibbs chains
  if(testCC){                        #retain random rows of W
    testCCn  <- 6
    testCCid <- sample(n,testCCn)
    wgibbs   <- array( data = 0, dim = c(ng, testCCn, S),
                       dimnames = list(NULL, NULL, snames) )
  }
  
  fgibbs <- matrix(0,ng,q1)
  colnames(fgibbs) <- fnames
  
  fbgibbs <- matrix(0,ng,q1*SO)
  colnames(fbgibbs) <- .multivarChainNames(fnames,snames[notOther])
  
  bgibbs <- matrix(0,ng,S*Q)
  colnames(bgibbs) <- .multivarChainNames(xnames,snames)
  
  covE <- cov( x%*%dCont )  # note that x is standardized
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)
  
  
  # unstandardized coefficients
  nsr <- length(standRows)
  
  tmp <- .getUnstandX(x, standRows, standMu[standRows,1],standMat[standRows,1],
                      intMat)
  S2U <- tmp$S2U
  xUnstand <- tmp$xu
  
  
  for(g in 1:ng){ ########################################################
 
    if(holdoutN > 0){
      tmp <- .getHoldLoHi( yh = yp[drop=F,holdoutIndex,], wh = w[drop=F,holdoutIndex,],
                           pl = plo[drop=F,holdoutIndex,],ph = phi[drop=F,holdoutIndex,],
                           eff = effMat[drop=F,holdoutIndex,], 
                           ymax = wmax, typeNames, cutg, ordCols)
      plo[holdoutIndex,] <- tmp$pl
      phi[holdoutIndex,] <- tmp$ph
    }
      
    tmp <- .param.fn(x, beta = bg[,notOther], Y = w[,notOther], otherpar)
    sg[notOther,notOther] <- tmp$sg
    bg[,notOther]         <- tmp$bg 
    otherpar              <- tmp$otherpar
    
    if(REDUCT){
      
      rndEff[,notOther]   <- tmp$rndEff
      sigmaerror          <- otherpar$sigmaerror
      kgibbs[g,notOther]  <- otherpar$K
      sgibbs[g,]          <- as.vector(otherpar$Z)
      sigErrGibbs[g]      <- otherpar$sigmaerror
      sg[sgOther]         <- otherpar$sigmaerror
    } else {
      sgibbs[g,] <- sg[Kindex]
    }
    
    alpha <- .sqrtRootMatrix(bg,sg,DIVIDE=T)
    
    if( 'OC' %in% typeCode ){
      tg   <- .gjamUpdateTheta(w,tg,cutLo,cutHi,ordCols,
                               holdoutN,holdoutIndex,minOrd,maxOrd) # var scale
      cutg <- .gjamCuts2theta(tg,ss = sg[ordCols,ordCols]) # corr scale
      breakMat[ordCols,1:lastOrd] <- cutg
      cgibbs[g,] <- as.vector( cutg[,-c(1,2,ncut)] )
      
      plo[,ordCols] <- cutg[cutLo]
      phi[,ordCols] <- cutg[cutHi]
      
   #   if(holdoutN > 0){
        
    #    ploHold <- plo[drop=F,holdoutIndex,]   --this will remain constant
    #    phiHold <- phi[drop=F,holdoutIndex,]
        
    #    tgHold <- .gjamUpdateTheta(wIn,tgHold,cutLo,cutHi,ordCols,
    #                           holdoutN=NULL,holdoutIndex=NULL,minOrd,maxOrd)
    }
    
    muw   <- x%*%bg
    
    tmp   <- .updateW( x, w, y, muw, sg, alpha, cutg, plo, phi, 
                       rndEff, sigmaerror, wHold )
    w     <- tmp$w
    yp    <- tmp$yp
    plo   <- tmp$plo
    phi   <- tmp$phi
    wHold <- tmp$wHold    #values for w if not held out
    
    if(testCC)wgibbs[g,,] <- w[testCCid,]
    
    setTxtProgressBar(pbar,g)
    
    if(mmiss > 0)y[ymiss] <- yp[ymiss]
    
    sinv <- .invertSigma(sg[notOther,notOther],sigmaerror,otherpar,REDUCT)
    
    if(nmiss > 0){
      x[xmiss] <- .imputX_MVN(x,w[,notOther],bg[,notOther],xmiss,sinv,xprior=xprior,
                              xbound=xbound)[xmiss]
      tmp    <- .getUnstandX(x, standRows, standMu[standRows,1],
                             standMat[standRows,1], intMat)
      S2U    <- tmp$S2U
      XX     <- crossprod(x)
      IXX    <- solve(XX)
    }
    
    bgs <- bg                        # unstandardize for X
    if(length(standRows) > 0)bgs <- S2U%*%muw
    
    bgibbs[g,] <- bgs
    
    if(TRAITS){
 
      Atrait <- bgs%*%t(specTrait[,colnames(yp)])
      Strait <- specTrait[,colnames(yp)]%*%sg%*%t(specTrait[,colnames(yp)])
      agibbs[g,] <- Atrait
      mgibbs[g,] <- Strait
      
   #   if(length(missTrait) > 0){
   #     yw     <- sweep(yp,1,rowSums(yp),'/')
   #     yw[yw <= 0]   <- 0
   #     yw[is.na(yw)] <- 0
   #     Ttrait <- .gjamPredictTraits(yw,specTrait[,colnames(yp)], traitTypes)
   #     tmp <- t(Ttrait)%*%yw
   #   }
      
    }
    
    if( PREDICTX & length(predXcols) > 0 ){
      
      ww <- w
      if(holdoutN > 0) ww[holdoutIndex,] <- wHold
        
      if( length(isNonLinX) > 0 ){
        
        xpred <- .predictY2X_nonLinear(xpred, yy=ww[,notOther],
                                       bb=bg[,notOther],ss=sg[notOther,notOther],
                                       priorIV = priorXIV,priorX=priorX,
                                       predCols=isNonLinX,isInt,intMat,
                                       isFactor,factorList, contrast = contrast,
                                       lox, hix)$x
        }
      
      if( length(px) > 0 ){
        
        wn <- which(!is.finite(xpred),arr.ind=T)
        if(length(wn) > 0){
          tmp <- matrix(priorX,Q,nrow(wn))
          xpred[wn[,1],] <- t(tmp)
        }
        xpred[,px] <- .predictY2X_linear(xpred, yy=ww[,notOther], bb=bg[,notOther],
                                         ss=sg[notOther,notOther], sinv = sinv,
                                         priorIV = priorXIV, 
                                         priorX=priorX,predCols=px, 
                                         REDUCT=REDUCT, lox, hix)[,px]
        wn <- which(!is.finite(xpred),arr.ind=T)
        if(length(wn) > 0){
          tmp <- matrix(priorX,Q,nrow(wn))
          xpred[wn[,1],] <- t(tmp)
        }
      }
      
      if( length(linFactor) > 0 ){
        
        xtmp <- xpred

        # predict all factors
        xtmp[,findex] <- .predictY2X_linear(xpred, yy=ww[,notOther], 
                                            bb=bg[,notOther],
                                  ss=sg[notOther,notOther], sinv = sinv,
                                  priorIV = priorXIV, 
                                  priorX=priorX,predCols=findex, 
                                  REDUCT=REDUCT, lox, hix)[,findex]
        
        for(k in 1:length(linFactor)){
          
          mm  <- linFactor[[k]]
          tmp <- xtmp[,mm]
          
          tmp[,1] <- 0
          ix  <- apply(tmp,1,which.max)   
          
          tmp <- tmp*0
        
          tmp[ cbind(1:n,ix) ] <- 1
          tmp <- tmp[,-1,drop=F]
          xpred[,mm[-1]] <- tmp
        }
      }
      xpred[,1] <- 1
    }
    
    # Fmatrix, bg is mostly standardized by x, bgg is completely standardized
    
    bgg <- bg[,notOther]
    if(!is.null(notStandard))bgg[notStandard,] <- bgg[notStandard,]*
      standMat[notStandard,notOther]
    
    agg <- .sqrtRootMatrix(bgg,sg[notOther,notOther],DIVIDE=T)  #cor-stand scale
    
    if(nfact > 0){
      agg <- lCont%*%agg    #standardized for x and cor scale for y
      for(k in 1:nfact){
        fk  <- facList2[[k]]
        mua <- colMeans(agg[drop=F,fk,])
        nl  <- length(fk)
        agg[fk,] <- agg[fk,] - matrix(mua,nl,SO,byrow=T)
      }
    } else {
      agg <- agg[drop=F,-1,]
    }
    
    egg         <- lCont%*%bgg          #standardized for x, not cor for y
    fsens       <- egg%*%sinv%*%t(egg)
    fgibbs[g,]  <- sqrt(diag(fsens))
    fbgibbs[g,] <- agg
    
    if(g > burnin){
      
      ntot   <- ntot + 1
      ypred  <- ypred + yp
      ypred2 <- ypred2 + yp^2
      sumDev <- sumDev - 2*sum(.dMVN(w[,notOther],muw[,notOther],
                                     sg[notOther,notOther], log=T ) )
      yerror <- yerror + (yp - y)^2
      
      fmat <- fmat + fsens
      
      sMean  <- sMean + sg
      
      wpred  <- wpred + w
      wpred2 <- wpred2 + w^2
      
      if(RICHNESS){
        yy <- yp[,notRichness,drop=F]
        if(length(notRichness) == 0)yy <- yp
        yy[yy > 0] <- 1
        yy[yy <= 0] <- 0
        richFull <- .add2matrix(rowSums(yy),richFull)
        richness <- .add2matrix(rowSums(yy*wrich[,notRichness,drop=F]),
                                richness)  # only for non-missing
      }
      
      if(mmiss > 0){
        ymissPred[ymiss]  <- ymissPred[ymiss] + y[ymiss]
        ymissPred2[ymiss] <- ymissPred2[ymiss] + y[ymiss]^2
      }
      
      if(PREDICTX & length(predXcols) > 0){
        predx  <- predx + xpred
        predx2 <- predx2 + xpred^2
      }
      
      ess[notOther,notOther]  <- .cov2Cor( t(agg)%*%covE%*%agg ) 
      emat[notOther,notOther] <- emat[notOther,notOther] + ess[notOther,notOther]
      
      lo[ ess < 0 ] <- lo[ ess < 0 ] + 1
      hi[ ess > 0 ] <- hi[ ess > 0 ] + 1
      
      ess[notOther,notOther] <- ginv(ess[notOther,notOther])
      
      lm[ ess < 0 ] <- lm[ ess < 0 ] + 1  # neg values
      hm[ ess > 0 ] <- hm[ ess > 0 ] + 1  # pos values
      
      if(TRAITS){
        yw     <- sweep(yp,1,rowSums(yp),'/')
        yw[yw <= 0]   <- 0
        yw[is.na(yw)] <- 0
        Ttrait <- .gjamPredictTraits(yw,specTrait[,colnames(yp)], traitTypes)
        tpred  <- tpred + Ttrait
        tpred2 <- tpred2 + Ttrait^2
      }
    }
  }     ################# end gibbs loop ####################
   
  otherpar$S <- S 
  otherpar$Q <- Q
  otherpar$snames <- snames
  otherpar$xnames <- xnames
  
  if(RICHNESS){
    richNonMiss <- richness/ntot            #only non-missing plots
    yr  <- as.matrix(ydata[,notRichness]) 
    yr[yr > 0] <- 1
    yr <- rowSums(yr,na.rm=T)
    vv  <- matrix(as.numeric(colnames(richNonMiss)),n,ncol(richNonMiss),byrow=T)
    rmu <- rowSums( vv * richNonMiss )/rowSums(richNonMiss)
    rsd <- sqrt( rowSums( vv^2 * richNonMiss )/rowSums(richNonMiss) - rmu^2)
    
    vv  <- matrix(as.numeric(colnames(richFull)),n,ncol(richFull),byrow=T)
    rfull <- rowSums( vv * richFull )/rowSums(richFull)
    richness <- cbind(yr, rmu, rsd, rfull )
    colnames(richness) <- c('obs','predMu','predSd','predNotMissing')
  }
    
  if(mmiss > 0){
    ymissPred[ymiss]  <- ymissPred[ymiss]/ntot
    ymissPred2[ymiss] <- sqrt(ymissPred2[ymiss]/ntot - ymissPred[ymiss]^2)
  }
  
  xunstand    <- .getUnstandX(x, standRows, standMu[standRows,1],
                         standMat[standRows,1], intMat)$xu
  
  rmspeBySpec <- sqrt( colSums(yerror)/ntot/n)
  rmspeAll    <- sqrt( sum(yerror)/ntot/n/S )
  
  sMean <- sMean/ntot
  
  
  betaMu <- colMeans(bgibbs[burnin:ng,])    #unstandardized
  betaSe <- apply(bgibbs[burnin:ng,],2,sd)
  betaMu <- matrix(betaMu,Q,S)
  betaSe <- matrix(betaSe,Q,S)
  colnames(betaMu) <- colnames(betaSe) <- snames
  rownames(betaMu) <- rownames(betaSe) <- xnames
  
  bplus  <- which( apply(bgibbs,2,quantile,.025) > 0)
  bminus <- which( apply(bgibbs,2,quantile,.975) < 0)
  
  fBetaMu <- colMeans(fbgibbs[burnin:ng,])
  fBetaMu <- matrix(fBetaMu,q1,SO)
  fBetaSd <- apply(fbgibbs[burnin:ng,],2,sd)
  fBetaSd <- matrix(fBetaSd,q1,SO)
  
  fMu <- colMeans(fgibbs[burnin:ng,,drop=F])
  fSe <- apply(fgibbs,2,sd)
  
  rownames(fBetaMu) <- rownames(fBetaSd) <- names(fMu)
  colnames(fBetaMu) <- colnames(fBetaSd) <- snames[notOther]
  
  yMu <- ypred/ntot
  ySd <- sqrt(ypred2/ntot - yMu^2)
  cMu <- cuts
  cSe <- numeric(0)
  
  wMu <- wpred/ntot
  wpp <- pmax(0,wpred2/ntot - wMu^2)
  wSd <- sqrt(wpp)
  
  
  # note: on latent w scale
  meanDev <- sumDev/ntot
  betaS   <- solve(crossprod(x))%*%crossprod(x,xunstand)%*%betaMu   #standardized
  
  pd  <- meanDev - 2*sum(.dMVN(wMu[,notOther],x%*%betaS[,notOther],
                               sMean[notOther,notOther], log=T ) )
  DIC <- pd + meanDev
  
  yscore <- colSums( .getScoreNorm(y[,notOther],yMu[,notOther],
                                   ySd[,notOther]^2),na.rm=T )  # gaussian w
  
  
  chains <- list( sgibbs = sgibbs, bgibbs = bgibbs) 

  if(REDUCT) chains <- append(chains,list(kgibbs = kgibbs,
                                          sigErrGibbs = sigErrGibbs))
  if(testCC) chains <- append(chains,list(wgibbs = wgibbs))
  
  xscore <- xpredMu <- xpredSd <- NULL
  standX <- xmissMu <- xmissSd <- NULL
  if(PREDICTX){
    xpredMu <- predx/ntot
    xpredSd <- predx2/ntot - xpredMu^2
    xpredSd[xpredSd < 0] <- 0
    xpredSd <- sqrt(xpredSd)
    
    xpredMu <- .getUnstandX(xpredMu, standRows, standMu[standRows,1],
                            standMat[standRows,1], intMat)$xu
    xpredSd[,standRows] <- xpredSd[,standRows]*
      matrix( standMat[standRows,1], n, length(standRows),byrow=T ) 
    
    if(Q == 2)xscore <- mean( .getScoreNorm(x[,2],xpredMu[,2],xpredSd[,2]^2) )
    if(Q > 2) xscore <- colMeans( .getScoreNorm(x[,-1],xpredMu[,-1],xpredSd[,-1]^2) )
  }

  if(length(standRows) > 0){                #unstandardize
    if(length(xmiss) > 0){
      wmm <- which( xnames[xmiss[,2]] %in% standRows )
      if(length(wmm) > 0){
        xmissMu <- x[xmiss]
        xmissSd[xmiss[wmm,]] <- xmissSd[xmiss[wmm,]]*standMat[xmiss[wmm,2],1]
      }
    }
    standX <- cbind(standMu[,1],standMat[,1])
    colnames(standX) <- c('xmean','xsd')
    rownames(standX) <- rownames(standMat)
  }

  # betaSens, sigma and R
  
  ns <- 500
  simIndex <- sample(burnin:ng,ns,replace=T)
  
  tmp <- .expandSigmaChains(snames, sgibbs, otherpar, simIndex=simIndex,
                            sigErrGibbs, kgibbs, REDUCT)
  corMu <- tmp$rMu; corSe <- tmp$rSe
  sigMu  <- tmp$sMu; sigSe  <- tmp$sSe
  
  
  whichZero <- which(lo/ntot < ematAlpha & 
                       hi/ntot < ematAlpha,arr.ind=T) #not different from zero
  whConZero <- which(lm/ntot < ematAlpha & 
                       hm/ntot < ematAlpha,arr.ind=T)
  
  ematrix  <- emat/ntot
  fmatrix  <- fmat/ntot
  
  tMu <- tSd <- tMuOrd <- btMu <- btSe <- stMu <- stSe <- numeric(0)
  
  if(TRAITS){
    
    tMu <- tpred/ntot
    tSd <- sqrt(tpred2/ntot - tMu^2)
    wo  <- which(traitTypes == 'OC')    #predict ordinal scores
    M   <- ncol(tMu)
    
    if(length(wo) > 0){
      tMuOrd <- tMu*0
      for(j in wo)tMuOrd[,j] <- round(tMu[,j],0) - 1
      tMuOrd <- tMuOrd[,wo]
    }
    
    tmp <- .processPars(agibbs)$summary
    btMu <- matrix(tmp[,'estimate'], Q, M)
    btSe <- matrix(tmp[,'se'], Q, M)
    
    tmp <- .processPars(mgibbs)$summary
    stMu <- matrix(tmp[,'estimate'],M,M)
    stSe <- matrix(tmp[,'se'],M,M)
    
    rownames(btMu) <- rownames(btSe) <- colnames(x)
    colnames(btMu) <- colnames(btSe) <- rownames(stMu) <- colnames(stMu) <- 
      rownames(stSe) <- colnames(stSe) <- tnames
    
    chains <- append( chains,list('agibbs' = agibbs))
    chains <- append( chains, list('mgibbs' = mgibbs) ) 
  }
  if('OC' %in% typeNames){
    nk  <- length(ordCols)
    tmp <- .processPars(cgibbs)$summary
    cMu <- matrix(tmp[,'estimate'],nk,ncut-3)
    cSe <- matrix(tmp[,'se'],nk,ncut-3)
    colnames(cMu) <- colnames(cSe) <- cnames[-c(1,2,ncut)]
    rownames(cMu) <- rownames(cSe) <- snames[ordCols]
    breakMat[ordCols,c(3:(3+(ncol(cMu))-1))] <- cMu
    chains <- c(chains,list(cgibbs = cgibbs))
  }
  
  chains <- c(chains, list(fgibbs = fgibbs, fbgibbs = fbgibbs) )
  
  if('PA' %in% typeNames){
    zMu <- yMu
    zSd <- ySd
  }
  
  colnames(betaMu)   <- colnames(betaSe) <- snames
  rownames(betaMu)   <- rownames(betaSe) <- xnames
  
  
  parameterTables <- list(tMuOrd = tMuOrd,tMu = tMu,tSd = tSd,
                          fMu = fMu, fSe = fSe, cutMu = cMu, cutSe = cSe, 
                          betaMu = betaMu, betaSe = betaSe, betaPlus = bplus,
                          betaMinus = bminus, corMu = corMu, corSe = corSe, 
                          fBetaMu = fBetaMu, fBetaSd = fBetaSd,
                          sigMu = sigMu, sigSe = sigSe, ematrix = ematrix, 
                          fmatrix = fmatrix,
                          betaTraitMu = btMu, betaTraitSe = btSe, 
                          sigmaTraitMu = stMu, sigmaTraitSe = stSe)
  parameterTables <- parameterTables[ sort( names(parameterTables) )]
  
  modelSummary <- append(modelSummary,
                         list(typeNames = typeNames, DIC = DIC, yscore = yscore, 
                              xscore = xscore, rmspeAll = rmspeAll,
                              rmspeBySpec = rmspeBySpec,
                              standMat = standMat, standRows = standRows,
                              xpredMu = xpredMu, xpredSd = xpredSd,
                              ypredMu = yMu, ypredSd = ySd, wMu = wMu, wSd = wSd, 
                              xdataNames = xdataNames,isNonLinX = isNonLinX,
                              ydataNames = ydataNames, isInt = isInt,
                              ematAlpha = ematAlpha, richness = richness,
                              whichZero = whichZero, whConZero = whConZero))
  
  modelSummary <- modelSummary[ sort( names(modelSummary) )]
  
  
  all <- list(effort = effort, REDUCT = REDUCT, 
              otherpar = otherpar, xdata = xdata, parameterTables = parameterTables,
       modelList = modelList, missingIndex = xmiss, xmissMu = xmissMu, 
       xmissSd = xmissSd, chains = chains, x = xunstand, y = y, holdoutIndex = holdoutIndex,
       yMissMu = ymissPred, yMissSd = ymissPred2, ymiss = ymiss, 
       modelSummary = modelSummary, breakMat = breakMat, standX = standX,
       censor = censor, TRAITS = TRAITS, traitList = traitList)
  all$call <- match.call()
  all <- all[ sort(names(all)) ]
  
  class(all) <- ".gjam"
  
 # .print.gjam(all)
  
  all
}

summary..gjam <- function(object,...){
  
  beta   <- object$parameterTables$betaMu
  oo <- grep('other',colnames(beta))
  if(length(oo) > 0)beta <- beta[,-oo]
  rb <- rownames(beta)
  cb <- colnames(beta)
 
  Sensitivity <- object$parameterTables$fMu
  StdErr      <- object$parameterTables$fSe
  sens <- data.frame(Sensitivity, StdErr)
  
  RMSPE       <- object$modelSummary$rmspeBySpec
  if(length(oo) > 0)RMSPE <- RMSPE[-oo]
  bb   <- signif(rbind(beta,RMSPE),3)
  
  cc <- as.vector( signif(beta, 3) )
  ss <- object$parameterTables$betaSe
  if(length(oo) > 0)ss <- ss[,-oo]
  ss <- as.vector( signif(ss, 3) )
  rr <- as.vector( outer(rb,cb,paste,sep='_') )
  TAB <- data.frame(Estimate = cc, StdErr = ss)
  rownames(TAB) <- rr
  
  cat("\nSensitivity by predictor variables f:\n")
  print( signif(sens, 3) )
  
  cat("\nCoefficient matrix B:\n")
  print(bb)
  cat("\nRMSPE by individual responses (columns).\n")
  
  cat("\nCoefficient matrix B with SEs:\n")
  
  Signif <- rep('',nrow(TAB))
  bp <- object$parameterTables$betaPlus
  bm <- object$parameterTables$betaMinus
  bs <- object$parameterTables$betaMu*0
  bs[bp] <- 1
  bs[bm] <- -1
  if(length(oo) > 0){
    bs <- bs[,-oo]
  }
  bs <- as.vector(bs)
  
  Signif[bs == 1] <- '+'
  Signif[bs == -1] <- '-'
  
  TAB <- data.frame(TAB,Signif)
  
  print(TAB)
  cat("\nThird column indicates if 95% posterior distribution contains zero.\n")
  
  res <- list(DIC=object$modelSummary$DIC, Sensitivity = sens, 
              Coefficients=TAB, beta = beta)
  class(res) <- "summary.gjam"
  invisible(res) 
}

print..gjam <- function(x, ...){
  
  cat("Call:\n")
  print(x$call)
  
  cat("\nDIC:\n")
  print( round(x$modelSummary$DIC,0) )

  cat("\nCoefficients:\n")
  print( signif(x$parameterTables$betaMu, 3) )
  
  cat("\nStandard errors:\n")
  print( signif(x$parameterTables$betaSe, 3) )
}

.getSigTable <- function(chain, SS, QQ, xn, sn){
  
  bci  <- apply(chain,2,quantile,c(.025,.975))
  tmp  <- .between(rep(0,SS*QQ),bci[1,],bci[2,],OUT=T)
  ii <- rep(' ',SS*QQ)
  ii[tmp[bci[1,tmp] < 0]] <- '-'
  ii[tmp[bci[2,tmp] > 0]] <- '+'
  
  bTab <- data.frame( matrix(ii,QQ,SS) )
  colnames(bTab) <- sn
  rownames(bTab) <- xn
  
  bTab <- data.frame( t(bTab) )
  
  bTab
}

.getPlotLayout <- function(np){
  
  # np - no. plots
  
  if(np == 1)return( c(1,1) )
  if(np == 2)return( c(1,2) )
  if(np == 3)return( c(1,3) )
  if(np <= 4)return( c(2,2) )
  if(np <= 6)return( c(2,3) )
  if(np <= 9)return( c(3,3) )
  if(np <= 12)return( c(3,4) )
  if(np <= 16)return( c(4,4) )
  if(np <= 20)return( c(4,5) )
  if(np <= 25)return( c(5,5) )
  if(np <= 25)return( c(5,6) )
  return( c(6,6) )
}

.gjamPlot <- function(output, plotPars){
  
  sdScaleX  <- plotAllY     <- sdScaleY <- ylog <- TRAITS <- F
  GRIDPLOTS <- SAVEPLOTS    <- REDUCT   <- TV   <- F
  sigOnly   <- PREDICTX     <- betaGrid <- PLOTY <- PLOTX  <- T 
  specLabs  <- SMALLPLOTS   <- T
  omitSpec   <- trueValues  <- censor <- otherpar <- NULL
  traitList  <- specByTrait <- typeNames <- classBySpec <- 
    x <- y   <- burnin      <- richness <- betaTraitMu <-  
    corSpec  <- cutMu       <- ypredMu <- DIC <- yscore <- missingIndex <- 
    xpredMu  <- plotByTrait <- tMu <- tMuOrd <- traitTypes <- 
    isFactor <- betaMu      <- corMu <- modelSummary <- NULL
  unstandardX  <- NULL # vector of std devs of x, by name, to unstandardize
  ematAlpha     <- .95  # threshold prob that a covariance/edge is not zero
  ematrix      <- NULL
  eCont <- modelList <- NULL
  
  ncluster <- min(c(4,ncol(output$y)))
  
  corLines <- T
  cex <- 1
  holdoutIndex <- numeric(0)
  clusterIndex <- clusterOrder <- numeric(0)
  parameterTables <- NULL

  outfolder <- 'gjamOutput'
  outfile   <- character(0)

  width <- height <- 3

  oma <- c(1,1,0,0)
  mar <- c(1,1,1,0)
  tcl <- -0.1
  mgp <- c(0,0,0)

  specColor <- traitColor <- textCol <- 'black'
  
  for(k in 1:length(output))assign( names(output)[k], output[[k]] )
  for(k in 1:length(modelList))assign( names(modelList)[k], modelList[[k]] )
  for(k in 1:length(modelSummary))assign( names(modelSummary)[k], modelSummary[[k]] )
  for(k in 1:length(parameterTables))assign( names(parameterTables)[k], 
                                             parameterTables[[k]] )
  for(k in 1:length(plotPars))assign( names(plotPars)[k], plotPars[[k]] )
  if( !is.null(traitList) ){
    TRAITS <- T
    for(k in 1:length(traitList))assign( names(traitList)[k], traitList[[k]] )
  }
  if( 'trueValues' %in% names(plotPars) ){
    TV <- T
    for(k in 1:length(trueValues))assign( names(trueValues)[k], trueValues[[k]] )
    
    matchTrue <- match(colnames(betaMu),colnames(beta))
    beta      <- beta[,matchTrue]
    sigma     <- sigma[matchTrue,matchTrue]
    corSpec   <- corSpec[matchTrue,matchTrue]
  }
  
  if(length(xpredMu) == 0)PREDICTX <- F
  if(!PREDICTX)PLOTX <- F
  
  if(!SMALLPLOTS){
    oma <- c(0,0,0,0)
    mar <- c(4,4,2,1)
    tcl <- -0.5
    mgp <- c(3,1,0)
  }
  if(SAVEPLOTS){
    ff <- file.exists(outfolder)
    if(!ff)dir.create(outfolder)
  }
  
  chainNames <- names(output$chains)
  
  allTypes <- unique(typeNames)
  
  ntypes   <- length(allTypes)
  typeCode <- match(typeNames,allTypes)
  specs    <- rownames(classBySpec)
  Q        <- ncol(x)
  nhold   <- length(holdoutIndex)
  ncut    <- ncol(classBySpec) + 1
  S       <- ncol(y)
  n       <- nrow(y)
  snames  <- colnames(y)
  xnames  <- colnames(x)
  ng      <- nrow(chains$bgibbs)
  gindex  <- burnin:ng
  
  if(S > 10)corLines <- F
  if(S < 8){
    if(GRIDPLOTS)message('no GRIDPLOTS if S < 8')
    GRIDPLOTS <- F
  }
  
  if(length(specColor) == 1)specColor <- rep(specColor, S)
  boxCol    <- .getColor(specColor,.4)
  
  for(k in 1:length(chains)){         # remove burnin
    if(length(chains[[k]]) > 0){
      if(is.matrix(chains[[k]])){
        chains[[k]] <- chains[[k]][gindex,]
      } else {
        chains[[k]] <- matrix(chains[[k]][gindex])
      }
    }
  }
  
  other    <- grep('other',colnames(y))
  omit     <- c(which(colnames(y) %in% omitSpec),other)
  notOther <- notOmit <- c(1:S)
  if(length(other) > 0)notOther <- notOther[!notOther %in% other]
  SO       <- length(notOther)
  if(length(omit) > 0)notOmit <- notOmit[-omit]
  SM       <- length(notOmit)
  
  snames  <- colnames(y)
  xnames  <- colnames(x)
  
  cvec    <- c('black','brown','orange')
  if(ntypes > 4)cvec <- c(cvec,'green','blue')
  colF    <- colorRampPalette(cvec)
  
  ## richness prediction
  
  xSd <- sqrt( diag(cov(x)) )
  
  if( !TRAITS & !is.null(richness) ){
    
    xlimit <- range(richness[,1])
    
    if(diff(xlimit) > 0){
      
      if(SAVEPLOTS)pdf( file=.outFile(outfolder,'richness.pdf') )
      
      par(mfrow=c(1,1),bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
      
      npp    <- length(richness[,1])/20
      ylimit <- range(richness[,2])
      
      tmp <- .bins4data(richness[,1],nPerBin=npp,breaks=NULL,LOG=F)
      breaks <- tmp$breaks
      bins   <- tmp$bins
      nbin   <- tmp$nbin
      
      ncc    <- max( c(100,max(richness[,1])/20) )
      xy     <- .gjamBaselineHist(richness[,1],bins=bins,nclass=ncc)
      xy[2,] <- ylimit[1] + .3*xy[2,]*diff(ylimit)/max(xy[2,])
      plot(xy[1,],xy[2,],col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
           xlab=' ',ylab='Predicted')
      polygon(xy[1,],xy[2,],border='brown',col='tan')
      
      fill <- .getColor('blue',.2)
      
      .plotObsPred(obs=richness[,1],yMean=richness[,2],
                   nPerBin=npp, wide=.6, MEDIAN = F, fill=fill, 
                   box.col='darkblue', 
                   xlabel='Observed',ylabel='Predicted', POINTS=F,
                   add=T)
      
      abline(0,1,lty=2)
      abline(h=mean(richness[,1]),lty=2)
      .plotLabel(' Richness (no. present)',cex=1.2,above=T)
      
      if(!SAVEPLOTS){
        readline('no. species may not vary much -- return to continue ')
      } else {
        dev.off( )
      }
    } 
  }
  
  #######################################
  
  tmp <- .omitChainCol(chains$bgibbs,'other')
  omitBC <- tmp$omit
  keepBC <- tmp$keep
  
  ns <- min(output$ng - output$burnin,500)
  simIndex <- sample(nrow(chains$sgibbs),ns,replace=T)
  
  tmp <- .expandSigmaChains(snames, chains$sgibbs, otherpar, simIndex, 
                            chains$sigErrGibbs, chains$kgibbs, 
                            REDUCT)
  corMu <- tmp$rMu; corSe <- tmp$rSe; sigMu  <- tmp$sMu; sigSe  <- tmp$sSe
  
  if(REDUCT){
    sigmaerror <- mean(chains$sigErrGibbs)
    sinv <- .invertSigma(sigMu,sigmaerror,otherpar,REDUCT)
  } else {
    sinv <- solve(sigMu[notOther,notOther])
  }
  
  bgibbsShort    <- chains$bgibbs[simIndex,]
  sgibbsShort    <- tmp$chainList$schain      #lower.tri with diagonal
  rgibbsShort    <- tmp$chainList$cchain
 
  if(REDUCT){
    kgibbsShort  <- tmp$chainList$kchain
    otherpar       <- output$otherpar
  }
  
  SO <- length(notOther)
  
  fMat <- output$parameterTables$fmatrix
  
  betaLab   <- expression( paste('Coefficient matrix ',hat(bold(B))  ))
  corLab    <- expression( paste('Correlation matrix ',hat(bold(R))  ))
  cutLab    <- expression( paste('Partition matrix ',hat(bold(plain(P)))  ))
  
  AA <- F
  if(!SMALLPLOTS)AA <- T
  
  ################if true parameter values
  
  if(TV){
    
    mfcol <- c(1,2)
    if('OC' %in% typeNames)mfcol = c(2,2)
    
    if(SAVEPLOTS)pdf( file=.outFile(outfolder,'trueVsPars.pdf') )
    
    par(mfcol=mfcol,bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
    colF    <- colorRampPalette(c('darkblue','orange'))
    cols    <- colF(ntypes)
    
    if('beta' %in% names(trueValues)){
      cols <- colF(ntypes)
      if(length(beta) < 100){
        .gjamTrueVest(chains$bgibbs[,keepBC],true=beta[keepBC],
                      typeCode,allTypes,colors=cols,label = betaLab)
      } else {
        .plotObsPred(beta[,notOther],betaMu[,notOther],xlabel='true',
                     ylabel='estimate', nPerBin=length(beta)/10,
                     fill='lightblue',box.col=cols,POINTS=T,MEDIAN=F,add=F)
        abline(0,1,lty=2)
      }
    }
    
    if( 'corSpec' %in% names(trueValues) ){
      
      cols <- colF(2^ntypes)
      
      corTrue <- corSpec
      diag(corTrue) <- NA
      if(length(other) > 0){
        corTrue[other,] <- NA
        corTrue[,other] <- NA
      }
      
      cindex <- which(lower.tri(corSpec,diag=T))            #location on matrix
      pindex <- which(lower.tri(corSpec,diag=T),arr.ind=T)

      rindex <- which(is.finite(corTrue[cindex]))     #location in sgibbs
      cindex <- cindex[rindex]
      pindex <- pindex[rindex,]   
      
      cols <- colF(ntypes + ntypes*(ntypes-1)/2)
      
      rg <- rgibbsShort
      rg[rg == 1] <- NA
      
      xlim <- range(c(-.1,.1,corTrue[cindex]),na.rm=T)
      ylim <- range(c(-.1,.1,rg),na.rm=T)
      
      add <- F
      m   <- 1
      combNames <- character(0)
      combCols  <- numeric(0)
      
      box <- F
      
      for(k in 1:length(allTypes)){
        
        wk <- which(typeNames == allTypes[k])
        wk <- wk[wk %in% notOther]
        wp <- which(pindex[,1] %in% wk & pindex[,2] %in% wk)  
        
        if( length(wp) == 1 ){
          
          combNames <- c(combNames,allTypes[k])
          
          yci <- quantile( rgibbsShort[,rindex[wp]] ,c(.5,.025,.975))
          xmu <- corSpec[matrix(pindex[wp,],1)]
          if(!add){
            plot(xmu,yci[1],xlim=xlim,ylim=ylim,
                 pch=3,col=cols[m], xlab='true',ylab='')
            add <- T
          } else {
            points(xmu,yci[1],pch=3,col=cols[m])
          }
          lines( c(xmu,xmu),yci[2:3],col=cols[m],lwd=2)
        }
        
        if(length(wp) > 1){
          
          if(length(wp) < 100){
            .gjamTrueVest(rgibbsShort[,rindex[wp]],true=corSpec[cindex[wp]],
                          typeCode,allTypes,label=corLab,xlim=xlim,ylim=ylim,
                          colors=cols[m],legend=F,add=add)
          } else {
            box <- T
            .plotObsPred(corSpec[cindex[wp]],corMu[cindex[wp]],xlabel='true',
                         ylabel='estimate', fill='lightblue',
                         nPerBin=length(wp)/10,box.col=cols[m], POINTS=T,
                         MEDIAN=F,add=add)
            if(!add)abline(0,1,lty=2)
          }
          add <- T
          combNames <- c(combNames,allTypes[k])
          combCols  <- c(combCols,cols[m])
          m <- m + 1
        }
        
        if(k < length(allTypes)){
          
          for( j in (k+1):length(allTypes) ){
            
            wj <- which(typeNames == allTypes[j])
            wj <- wj[wj %in% notOther]
            wp <- which(pindex[,1] %in% wk & pindex[,2] %in% wj)
            
            if(length(wp) == 0){
              wp <- which(pindex[,2] %in% wk & pindex[,1] %in% wj)
            }
            
            if(length(wp) == 0)next
            
            if(length(wp) == 1){
              yci <- quantile( rgibbsShort[,rindex[wp]] ,c(.5,.025,.975))
              xmu <- corTrue[cindex[wp]]
              if(!add){
                plot(xmu,yci[1],xlim=xlim,ylim=ylim,
                     pch=3,col=cols[m])
              } else {
                points(xmu,yci[1],pch=3,col=cols[m])
              }
              lines( c(xmu,xmu),yci[2:3],col=cols[m],lwd=2)
              
            } else {
              if(!box){
                .gjamTrueVest(rgibbsShort[,rindex[wp]],
                             true=corTrue[cindex[wp]],
                             typeCode,allTypes,add=add,colors=cols[m],
                             legend=F, xlim=c(-.9,.9), ylim=c(-.9,.9))
              } else {
                .plotObsPred(corSpec[cindex[wp]],corTrue[cindex[wp]],
                             nPerBin=length(wp)/10,
                             box.col=cols[m], fill='white',POINTS=T,
                             MEDIAN=F,add=add)
              }
            }
            m <- m + 1
            
            mnames    <- paste(allTypes[k],allTypes[j],sep='-')

            combNames <- c(combNames,mnames)
            combCols  <- c(combCols,rep(cols[m],length(mnames)))
            add <- T
          }
          
        }
      }
      legend('topleft',combNames,text.col=cols,bty='n',ncol=3,cex=.7)
      
    }
    
    if('OC' %in% allTypes & 'cuts' %in% names(trueValues)){
      wc       <- c(1:ncol( cutMu )) + 2
      ctrue    <- cuts[,wc]
      wf       <- which(is.finite(ctrue*cutMu))
      cutTable <- .gjamTrueVest(chains$cgibbs[,wf],true=ctrue[wf],
                                typeCode,allTypes,colors='black',
                                label=cutLab,legend=F, add=F)
    }
    
    if(!SAVEPLOTS){
      readline('simulated beta, corSpec vs betaMu, corMu (95%) -- return to continue')
    } else {
      dev.off()
    }
  }
  
  ##################### partition for ordinal
  
  if('OC' %in% typeNames){
    
    if(SAVEPLOTS)pdf( file=.outFile(outfolder,'partition.pdf') )
    
    par( mfrow=c(1,1), bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp )
    
    wk <- which(typeNames == 'OC')
    nk <- length(wk)
    
    cgibbs <- chains$cgibbs
    
    cgibbs[!is.finite(cgibbs)] <- NA
    cc <- colSums(abs(cgibbs),na.rm=T)
    cg <- cgibbs[,cc > 0]
    
    # vnames <- rownames(cutMu)
    if('cuts' %in% names(trueValues))rownames(cuts) <- rownames(cutMu)
    
    c1 <- names(cc)[cc > 0]
    vnames <- sort(unique(matrix(unlist(strsplit(c1,'_')),ncol=2,byrow=T)[,1]))
    
    colc <- colF(ncol(cutMu))
    
    nk <- length(vnames)
    plot(0,0,xlim=c(0,max(cg,na.rm=T)),ylim=c(1,1+nk),cex=.1,
         xlab='Unit variance scale',
         ylab='Species',yaxt='n')
    .yaxisHorizLabs(vnames,at=c(1:nk))
    
    for(k in 1:length(vnames)){
      
      x1   <- 0
      ym   <- .5
      
      tmp <- .chains2density(cg,varName=vnames[k], cut=2.5)
      xt  <- tmp$x
      yt  <- tmp$y
      yt  <- .5*yt/max(yt)
      
      yt <- yt + k
      
      for(j in 1:nrow(xt)){
        
        if('cuts' %in% names(trueValues)){
          lines( rep(cuts[vnames[k],j+2],2),c(k,k+1),lty=2,col=colc[j],lwd=3)
        }
        
        xj <- c(xt[j,],xt[j,ncol(xt)],xt[j,1])
        yj <- c(yt[j,],k,k)
        
        x2 <- which.max(yj)
        xm <- .2*x1 + .8*xj[x2]
        
        polygon(xj,yj,border=colc[j],col=colc[j],lwd=2)
        if(k == length(vnames)) text(xm,ym+k,j,col=colc[j])
        x1 <- xj[x2]
      }
    }
    .plotLabel('Partition by species',above=T)
    
    if(!SAVEPLOTS){
      readline('cuts vs cutMu -- return to continue')
    } else {
      dev.off()
    }
  }  
  
  ############################
  
  rmspeAll <- sqrt( mean( (y[,notOther] - ypredMu[,notOther])^2 ) )
  
  eBySpec <- sqrt( colMeans( (y[,notOther]/rowSums(y[,notOther]) - 
                                ypredMu[,notOther]/rowSums(ypredMu[,notOther]))^2 ) )
  ordFit  <- order(eBySpec)
  
  score <- mean(yscore)
  
  fit <- signif( c(DIC,score,rmspeAll), 5)
  names(fit) <- c('DIC','score','rmspe')
  
  ################## predict y
  
  if(PLOTY){
    
    if(SAVEPLOTS)pdf( file=.outFile(outfolder,'yPred.pdf') )
    
    npp <- 0
    
    for(k in 1:length(allTypes)){
      wk    <- which(typeCode == k)
      if( length(censor) > 0 ){
        ncc <- 0
        if( typeNames[wk[1]] %in% names(censor) ){
          wm   <- which(names(censor) == typeNames[wk[1]])
          #    wall <- wm
          wnot <- wk
          for(m in wm){
            wnot <- wnot[!wnot %in% censor[[m]]$columns]
            npp  <- npp + 1
          }
          if(length(wnot) > 0)npp <- npp + 1
        } else {
          ncc <- ncc + 1
        }
      } else {
        ncc <- 1
      }
      npp <- npp + ncc
    }  
    
    mfrow <- .getPlotLayout(npp)
    par( mfrow=mfrow, bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp )
    
    ylab <- ' '
    mk   <- 0
    
    for(k in 1:length(allTypes)){
      
      wk    <- which(typeCode == k)
      wk    <- wk[wk %in% notOther]
      wkm   <- wk
      nk    <- nkm <- length(wk)
      censm <- NULL
      wm    <- wall <- 1
      CENS  <- F
      add   <- F
      
      if( length(censor) > 0 ){
        if( typeNames[wk[1]] %in% names(censor) ){
          CENS <- T
          wm   <- which(names(censor) == typeNames[wk[1]])
          wall <- wm
          wnot <- wk
          for(m in wm){
            wnot <- wnot[!wnot %in% censor[[m]]$columns]
          }
          if(length(wnot) > 0)wall <- c(wall,max(wall) + 1)
        }
      }
      
      for(m in wall){
        
        if(CENS){
          if(m %in% wm){
            censm <- censor[[m]]
            wkm   <- censor[[m]]$columns
          } else {
            censm <- NULL
            wkm   <- wnot
          }
          nkm <- length(wkm)
        }
        
        mk <- mk + 1
        
        y1 <- y[,wkm,drop=F]
        yp <- ypredMu[,wkm,drop=F]
        
        tmp <- .gjamPlotPars(type=typeNames[wk[1]],y1,yp,censm)
        y1 <- tmp$y1; yp <- tmp$yp; nbin <- tmp$nbin; nPerBin <- tmp$nPerBin
        vlines <- tmp$vlines; xlimit <- tmp$xlimit; ylimit <- tmp$ylimit
        breaks <- tmp$breaks; wide <- tmp$wide; LOG <- tmp$LOG; POINTS <- F
        MEDIAN <- tmp$MEDIAN
 
        log <- ''
        if(LOG)log <- 'xy'
        
        if(LOG){
          wn <- which(y1 > 0 & yp > 0)
          y1 <- y1[wn]
          yp <- yp[wn]
        }
        
        tmp <- .bins4data(y1,nPerBin=nPerBin,breaks=breaks,LOG=LOG)
        breaks <- tmp$breaks
        bins   <- tmp$bins
        nbin   <- tmp$nbin
        
        if(length(bins) > 0){
          breaks <- bins
          nPerBin <- NULL
        }
        
        if( !typeNames[wk[1]] %in% c('PA','CAT') ){
          ncc   <- max( c(100,max(y1)/20) )
          xy <- .gjamBaselineHist(y1,bins=bins,nclass=ncc)
          xy[2,] <- ylimit[1] + .3*xy[2,]*diff(ylimit)/max(xy[2,])
          xy[1,xy[1,] < xlimit[1]] <- xlimit[1]
          xy[2,xy[2,] < ylimit[1]] <- ylimit[1]
          
          plot(xy[1,],xy[2,],col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
               xlab='Observed',ylab='Predicted', log=log)
          polygon(xy[1,],xy[2,],border='tan',col='wheat')
          
        } else {
          y11 <- mean(y1)
          y00 <- 1 - y11
          x11 <- c(-.07,-.07,.07,.07,.93,.93,1.07,1.07,-.07)
          y11 <- c(0,y00,y00,0,0,y11,y11,0,0)
          plot(x11,y11,col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
               xlab=' ',ylab=ylab)
          polygon(x11,y11,border='tan',col='wheat')
        }
        abline(0,1,lty=2,lwd=3,col='brown')
        abline(h = mean(y1),lty=2,lwd=3,col='tan')
        
        add <- T
        
        if(nhold > 0){
          y1h <- y[holdoutIndex,wkm,drop=F]
          yph <- ypredMu[holdoutIndex,wkm,drop=F]
          points(y1h,yph,col='brown',
                 pch=21, bg='green',cex=.3)
        } 
        
        fill <- .getColor('blue',.3)
        
        if(xlimit[2] < max(bins))xlimit[2] <- max(bins) + 1
        
        tmp <- .plotObsPred(y1,yp,xlabel='Observed',ylabel=ylab,nbin=nbin, 
                            nPerBin=NULL, xlimit=xlimit,ylimit=ylimit,
                            breaks=bins, wide=wide, LOG=LOG, fill=fill, 
                            box.col='darkblue',POINTS=F, MEDIAN=MEDIAN, add=add)
        
        if(length(vlines) > 0)abline(v=vlines,lty=2)
        
        tf <- .gjamGetTypes(typeNames[wk[1]])$labels
        tf <- paste(letters[mk],tf, sep=') ')
        
        .plotLabel(tf,'topleft',above=AA)

      }
    }
    
    if(!SAVEPLOTS){
      readline('obs y vs predicted y -- return to continue ')
    } else {
      dev.off()
    }  ##########################
    
    # mean/variance in y
    
    PLOTEFFORT <- F
    
    if('CC' %in% typeNames & PLOTEFFORT){
      
      if(SAVEPLOTS)pdf( file=.outFile(outfolder,'effort.pdf') )
      
      par(mfrow=c(1,1),bty='n')
      ypredMu <- output$modelSummary$ypredMu
      ySd <- output$modelSummary$ypredSd
      
      n <- nrow(y)
      effort <- matrix(rowSums(y),n,S)
      cvv <- ySd/ypredMu
      www <- which(is.finite(effort) & is.finite(cvv) & cvv > 0)
      
      yl <- log10(quantile(cvv[www],c(.05,.999)))
      yl[1] <- floor(yl[1])
      yl[2] <- ceiling(yl[2])
      xl <- log10(range(effort[www]))
      xl[1] <- floor(xl[1])
      xl[2] <- ceiling(xl[2])
      
      .plotObsPred(effort[www],cvv[www],xlabel='Effort',ylabel='cv(Y)',
                   xlimit = 10^xl, ylimit = 10^yl,
                   breaks=10^seq(xl[1],xl[2],length=15), fill='lightblue', LOG=T,
                   box.col='darkblue',POINTS=F, add=F)
      
      if(!SAVEPLOTS){
        readline('mean/var in CC data -- return to continue ')
      } else {
        dev.off()
      }
    }
  }  #################
  
  nfact <- length(modelSummary$factorList)
  
  if( PLOTX & PREDICTX & length(output$modelSummary$xpredMu) > 0){
    
    noX <- character(0)
    colorGrad   <- colorRampPalette(c('white','brown','black'))
    
    if(nfact > 0){
      
      nn <- length(unlist(modelSummary$factorList)) # + nfact
      mmat <- matrix(0,nn,nn)
      mnames <- rep('bogus',nn)
      samples <- rep(0,nn)
      
      ib <- 1
      
      par(mfrow=c(1,1),bty='n')
      
      mm <- max(nfact,2)
      useCols <- colorRampPalette(c('brown','orange','darkblue'))(mm)
      textCol <- character(0)
      
      
      for(kk in 1:nfact){
        
        gname <- names( modelSummary$factorList )[[kk]]
        fnames <- modelSummary$factorList[[kk]]
        nx     <- length(fnames)
        if(nx < 1)next
        
        ie <- ib + nx - 1
        
        noX <- c(noX,fnames)
        
        cont <- contrasts[[kk]]
        
        refClass <- names (which( rowSums( cont ) == 0) )
        
        hnames <- matrix( unlist( strsplit(fnames,gname) ),nx,2,byrow=T)[,2]
     #   hnames <- c(refClass,hnames)
        knames <- c(paste(gname,'Ref',sep=''),fnames)
        xtrue     <- x[,fnames,drop=F]
     #   reference <- 1 - rowSums(xtrue)
     #   xtrue <- cbind(reference,xtrue)
        nx    <- ncol(xtrue)
     #   colnames(xtrue)[1] <- refClass
        
        xpred <- xpredMu[,fnames,drop=F]
   #     reference <- 1 - rowSums(xpred)
   #     xpred <- cbind(reference,xpred)
        cmat  <- matrix(0,nx,nx)
        colnames(cmat) <- hnames
        rownames(cmat) <- rev(hnames)
        #    wt <- apply(xtrue,1,which.max)
        for(j in 1:nx){
          wj <- which(xtrue[,j] == 1)
          cmat[,j] <- rev( colSums(xpred[drop=F,wj,])/length(wj) )
        }
        
        nb <- nn - ib + 1
        ne <- nn - ie + 1
        samples[ib:ie] <- colSums(xtrue)/n
        
        mmat[ne:nb,ib:ie] <- cmat
        mnames[ib:ie] <- hnames 
        
        textCol <- c(textCol,rep(useCols[kk],nx))
        
        ib <- ie + 1
      }
      
      colnames(mmat) <- mnames
      rownames(mmat) <- rev(mnames)
      
      names(textCol) <- rev( rownames(mmat) )
      
    #  mmat <- rbind(mmat,samples*0,samples)
    #  textCol <- c(textCol,'black')
        
      graphics.off()
      if(SAVEPLOTS)pdf( file=.outFile(outfolder,'xPredFactors.pdf' ) )
      par(mfrow=c(1,1),bty='n')
      
        .corPlot(mmat,slim=c(0,1),plotScale=.8, textCol = textCol,
                 PDIAG=F,corLines=T, tri='both',
                 specLabs = T, colorGrad = colorGrad,
                 textSize=1, new=F)
        if(nx > 2){
          mloc <- par('usr')
          text(mean(mloc[1:2]),mloc[3] + .03*diff(mloc[3:4]),'Observed')
          mtext('Predicted',side=4)
        }
        
        if(!SAVEPLOTS){
          readline('x inverse prediction, factors -- return to continue ')
        } else {
          dev.off()
        }
      }

    noplot <- c(1,grep(':',xnames),grep('^2',xnames,fixed=T))
    vnames <- xnames[-noplot]
    vnames <- vnames[!vnames %in% noX]
    
    if(length(vnames) > 0){
      
      if(SAVEPLOTS)pdf( file=.outFile(outfolder,'xPred.pdf') )
      
      ylab <- 'Predicted'
      xlab <- ' '
      mfrow <- .getPlotLayout( length(vnames) )
      par( mfrow=mfrow, bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp )
      
      missX <- missingIndex
      xmaxX <- apply(x,2,max,na.rm=T)
      
      k <- 0
      b <- 0
      
      for(j in 2:Q){
        
        if(!xnames[j] %in% vnames)next
        xlab <- 'Observed'
        
        k <- k + 1
        b <- b + 1
        if(k > mfrow[1])ylab <- ' '
        if(b == mfrow[2]){
          b <- 0
          xlab <- 'Observed'
        }
        
        x1 <- x[,j]
        x2 <- xpredMu[,j]
        
        type <- 'CON'
        for(kk in 1:length(modelSummary$factorList)){
          if( xnames[j] %in% modelSummary$factorList[[kk]] )type <- 'PA'
        }
        
        tmp <- .gjamPlotPars(type=type,x1,x2)
        y1 <- tmp$y1; yp <- tmp$yp; nbin <- tmp$nbin; nPerBin <- tmp$nPerBin
        vlines <- tmp$vlines; xlimit <- tmp$xlimit; ylimit <- tmp$ylimit
        breaks <- tmp$breaks; wide <- tmp$wide; LOG <- tmp$LOG; POINTS <- F
        MEDIAN <- tmp$MEDIAN
        
        LOG <- F
        
        if(nhold > 0){
          x1 <- x1[-holdoutIndex]
          x2 <- x2[-holdoutIndex]
          y1 <- y1[-holdoutIndex,,drop=F]
          yp <- yp[-holdoutIndex,,drop=F]
        }
        
        log <- ''
        if(LOG)log <- 'xy'
        
        if(LOG){
          wn <- which(y1 > 0 & yp > 0)
          y1 <- y1[wn]
          yp <- yp[wn]
        }
        
        tmp <- .bins4data(y1,nPerBin=nPerBin,breaks=breaks,LOG=LOG)
        breaks <- tmp$breaks
        bins   <- tmp$bins
        nbin   <- tmp$nbin
        
        if(length(bins) > 0){
          breaks <- bins
          nPerBin <- NULL
        }
        
        ncc   <- max( c(100,max(y1)/20) )
        xy <- .gjamBaselineHist(y1,bins=bins,nclass=ncc)
        xy[2,] <- ylimit[1] + .3*xy[2,]*diff(ylimit)/max(xy[2,])
        plot(xy[1,],xy[2,],col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
             xlab=' ',ylab=ylab)
        polygon(xy[1,],xy[2,],border='tan',col='wheat')
        
        
        abline(0,1,lty=2,lwd=3,col='brown')
        abline(h = mean(y1),lty=2,lwd=3,col='tan')
        
        add <- T
        
        if(nhold > 0){
          points(x[holdoutIndex,j],xpredMu[holdoutIndex,j],col='brown',
                 pch=21, bg='blue',cex=.4)
        } 
        
        
        fill <- .getColor('blue',.3)
        
        
        tmp <- .plotObsPred(y1,yp,xlabel='Observed',ylabel=ylab,nbin=nbin, 
                            nPerBin=nPerBin, xlimit=xlimit,ylimit=ylimit,
                            breaks=breaks, wide=wide, LOG=LOG, fill=fill, 
                            box.col='darkblue',POINTS=F, MEDIAN=MEDIAN, add=add)
        
        if(nhold > 0)points(x[holdoutIndex,j],xpredMu[holdoutIndex,j],
                            col='brown',cex=.3)
        
        if(length(missX) > 0){
          ww <- which(missX[,2] == j)
          if(length(ww) > 0){
            wz <- missX[ww,]
            if(!is.matrix(wz))wz <- matrix(wz,1)
            points(jitter(ww*0+xmaxX[j]),xpredMu[wz],cex=.6,col='blue')
          }
        }
        
        .plotLabel(paste(letters[j-1],xnames[j],sep=') '), above=AA)
      }
      
      if(!SAVEPLOTS){
        readline('x inverse prediction, covariates -- return to continue ')
      } else {
        dev.off()
      }
    }
  }   
  ######################
  
  if(plotAllY){
    
    if(SAVEPLOTS)pdf( file=.outFile(outfolder,'yPredAll.pdf') )
    
    np <- S   <- ncol(y)
    o   <- 1:S
    if(S > 16){
      np <- 16
      o <- sample(S)
    }
    
    mfrow <- .getPlotLayout(np)
    par(mfrow=mfrow, bty='n', oma=oma, mar=c(1,2,3,1), tcl= tcl, mgp=mgp)
    
    k <- 0
    add <- F
    
    for(j in o){
      
      censm <- NULL
      if( length(censor) > 0 ){
        if( typeNames[j] %in% names(censor) ){
          wjc <- which(names(censor) == typeNames[j])
          if(j %in% censor[[wjc]]$columns)censm <- censor[[wjc]]
        }
      }
      
      y1 <- y[,j]
      if(min(y1) == max(y1))next
      y2 <- ypredMu[,j]
      
      k <- k + 1
      if(k > 16)break
      
      tmp <- .gjamPlotPars(type=typeNames[j],y1,y2,censm)
      y1 <- tmp$y1; yp <- tmp$yp; nbin <- tmp$nbin; nPerBin <- tmp$nPerBin
      vlines <- tmp$vlines; xlimit <- tmp$xlimit; ylimit <- tmp$ylimit
      breaks <- tmp$breaks; wide <- tmp$wide; LOG <- tmp$LOG; POINTS <- F
      MEDIAN <- tmp$MEDIAN
      
      log <- ''
      if(LOG)log <- 'xy'
      
      if(LOG){
        wn <- which(y1 > 0 & yp > 0)
        y1 <- y1[wn]
        yp <- yp[wn]
      }
      
      LOG <- F
      
      tmp <- .bins4data(y1,nPerBin=nPerBin,breaks=breaks,LOG=LOG)
      breaks <- tmp$breaks
      bins   <- tmp$bins
      nbin   <- tmp$nbin
      
      
      if(length(bins) > 0){
        breaks <- bins
        nPerBin <- NULL
      }
      
      
      if( !typeNames[wk[1]] %in% c('PA','CAT') ){
        ncc   <- max( c(100,max(y1)/20) )
        xy <- .gjamBaselineHist(y1,bins=bins,nclass=ncc)
        xy[2,] <- ylimit[1] + .8*xy[2,]*diff(ylimit)/max(xy[2,])
        plot(xy[1,],xy[2,],col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
             xlab=' ',ylab=ylab)
        polygon(xy[1,],xy[2,],border='tan',col='wheat')
        
      } else {
        y11 <- mean(y1)
        y00 <- 1 - y11
        x11 <- c(-.07,-.07,.07,.07,.93,.93,1.07,1.07,-.07)
        y11 <- c(0,y00,y00,0,0,y11,y11,0,0)
        plot(x11,y11,col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
             xlab=' ',ylab=ylab)
        polygon(x11,y11,border='tan',col='wheat')
      }
      abline(0,1,lty=2,lwd=3,col='brown')
      abline(h = mean(y1),lty=2,lwd=3,col='tan')
      
      add <- T
      
      if(nhold > 0){
        points(y1[holdoutIndex],yp[holdoutIndex],col='brown',
               pch=21, bg='blue',cex=.4)
      } 
      
      fill <- .getColor('blue',.3)
      
      tmp <- .plotObsPred(y1,yp,xlabel='Observed',ylabel=ylab,nbin=nbin, 
                          nPerBin=nPerBin, xlimit=xlimit,ylimit=ylimit,
                          breaks=breaks, wide=wide, LOG=LOG, fill=fill, 
                          box.col='darkblue',POINTS=F, MEDIAN=MEDIAN, add=add)
      if(length(vlines) > 0)abline(v=vlines,lty=2)
      
      lab <- paste(letters[k],') ',colnames(y)[k],' - ', 
                   typeNames[k], sep='')
      
      .plotLabel( lab,above=T )
      abline(0,1,lty=2)
      abline(h = mean(y2),lty=2)
    }
    
    if(!SAVEPLOTS){
      readline('y prediction -- return to continue ')
    } else {
      dev.off()
    }
  }
  
  ##############traits
  
  if(TRAITS){
    
    if(SAVEPLOTS)pdf( file=.outFile(outfolder,'traitPred.pdf') ) # start plot
    
    tt <- grep('other',colnames(plotByTrait))
    if(length(tt) > 0)colnames(plotByTrait)[tt] <- colnames(specByTrait)[tt]
    
    print(colnames(plotByTrait))
               
    yy <- plotByTrait
    o  <- 1:ncol(yy)
    
    if(ncol(yy) > 16){
      
      rmspe <- sqrt( colSums( (plotByTrait - tMu)^2 )/n )
      o <- order(rmspe)[1:9]
      yy <- plotByTrait[,o]
    }
    
    mfrow <- .getPlotLayout(length(o))
    par(mfrow=mfrow, bty='n', oma=oma, mar=c(3,3,1,1), tcl= tcl, mgp=mgp)
    k <- 0
    
    for(j in o){
      
      add   <- F
      jname <- colnames(tMu)[j]
      
      k <- k + 1
      
      td <- plotByTrait[,jname]
      
      tjj   <- tMu[,j]
      wj <- which(colnames(tMuOrd) == jname)
      
      tmp <- .gjamPlotPars(type=traitTypes[j],td,tjj)
      y1 <- tmp$y1; yp <- tmp$yp; nbin <- tmp$nbin; nPerBin <- tmp$nPerBin
      vlines <- tmp$vlines; xlimit <- tmp$xlimit; ylimit <- tmp$ylimit
      breaks <- tmp$breaks; wide <- tmp$wide; LOG <- tmp$LOG; POINTS <- F
      MEDIAN <- tmp$MEDIAN
      
      if(nhold > 0){
        add <- T
        log <- ''
        if(LOG)log <- 'xy'
        plot(td[holdoutIndex],tjj[holdoutIndex],xlab=' ',ylab=ylab,
             xlim=xlimit,ylim=ylimit,col='grey',pch=21,bg='brown',cex=.4,log=log)
      } 
      
      tmp <- .plotObsPred(td,tjj,xlabel=' ',ylabel=ylab,nbin=nbin, 
                          nPerBin=nPerBin,
                          xlimit=xlimit,ylimit=ylimit,breaks=breaks,
                          wide=wide,LOG=LOG,
                          fill='grey',
                          POINTS=F,MEDIAN=MEDIAN,add=add)
      if(length(vlines) > 0)abline(v=vlines,lty=2)
      abline(0,1,lty=2)
      abline(h=mean(td,na.rm=T),lty=2)
      
      .plotLabel( paste(letters[k],') ',.traitLabel(jname),sep=''),above=AA )
    }
    
    if(!SAVEPLOTS){
      readline('predictive trait distributions -- return to continue ')
    } else {
      dev.off()
    }
  }
  
  ##############sensitivity 
  
  fgibbs <- chains$fgibbs
  if(!is.matrix(fgibbs)){
    fgibbs <- matrix(fgibbs)
    colnames(fgibbs) <- xnames[-1]
  }
  
  wc <- c(1:ncol(fgibbs))
  wx <- grep(':',colnames(fgibbs))
  wx <- c(wx, grep('^2',colnames(fgibbs), fixed=T) )
  if(length(wx) > 0)wc <- wc[-wx]
  
  wx <- grep('intercept',colnames(fgibbs))
  if(length(wx) > 0)wc <- wc[-wx]
  
  if(SAVEPLOTS)pdf( file=.outFile(outfolder,'sensitivity.pdf') ) # start plot
  
  xx   <- fgibbs[,wc,drop=F]
  tcol <- rep('black',ncol(xx))
  names(tcol) <- colnames(xx)
  
  if(nfact > 0){
    for(i in 1:nfact){
      im <- which(colnames(xx) %in% rownames(modelSummary$contrasts[[i]]))
      tcol[im] <- useCols[i]
    }
  }
  
  
#  tcol[ names(textCol) ] <- textCol
    
  par(mfrow=c(1,1),bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
  
  ord  <- order( colMeans(xx) )
  ylim <- c(min(xx),2*quantile(xx,.95))
  tmp <- .boxplotQuant( xx[,ord, drop=F], xaxt='n',outline=F, 
                        border=tcol[ord],whiskcol=tcol[ord],
                        boxfill=.getColor(tcol[ord],.4), 
                        pars = list(boxwex = 0.5, ylim=ylim), lty=1, log='y')
  mtext('Predictors in X',side=1,line=1)
  abline(h=0,lwd=2,col='grey')
  
  dy <- .05*diff(par()$yaxp[1:2])
  text(1:length(wc), dy + tmp$stats[5,],tmp$names,srt=90,pos=4,col=tcol[ord])
  
  sensLab   <- expression( paste('Sensitivity ',hat(bold(F))  ))
  
  .plotLabel(sensLab,above=AA, cex=.7)   
  
  if(!SAVEPLOTS){
    readline('sensitivity over full model -- return to continue ')
  } else {
    dev.off()
  }
  
  ######################  coefficient summary tables ############
  
  bTab   <- .getSigTable(chains$bgibbs,S, Q, xnames, snames) 
  
  q1    <- nrow(output$modelSummary$eCont)
  fnames <- rownames(output$modelSummary$eCont)
  bfTab <- .getSigTable(chains$fbgibbs,SO, q1, fnames, 
                        colnames(output$parameterTables$fBetaMu)) 
  
  bfCoeffTable <- .processPars(chains$fbgibbs,sigOnly=T)$summary
  sigFbeta     <- rownames(bfCoeffTable)
  
  bfSig <- chains$fbgibbs[,sigFbeta]
  
  
  #  bCoeffTable <- .processPars(chains$bgibbs[,keepBC],sigOnly=sigOnly)$summary
  #  sigBeta     <- rownames(bCoeffTable)
  
  bCoeffTable <- .processPars(chains$bgibbs[,keepBC],sigOnly=T)$summary
  sigBeta     <- rownames(bCoeffTable)
  bCoeffTable <- .processPars(chains$bgibbs[,keepBC],sigOnly=F)$summary
  
  if(length(sigBeta) == 0)sigBeta <- c(1:ncol(chains$bgibbs))
  
  scaleNote <- 'W/X scale'
  
  betaSig <- chains$bgibbs[,sigBeta]
  
  summaryCoeffs <- list(betaSig = bTab, fBetaSig = bfTab, 
                        betaCoeff = bCoeffTable, fBetaCoeff = bfCoeffTable)
  
#  if(sdScaleY){
#    ssi <- 1:(S*(S+1)/2)
 #   svi <- matrix( 0,S,S)
#    svi[ lower.tri(svi, diag=T) ] <- ssi
#    svi <- diag(svi)
#    svi <- rep(svi,each=Q)
#    bgibbsShortCor <- bgibbsShort/sqrt( sgibbsShort[,svi] )
#    betaSig <- bgibbsShortCor[,sigBeta]
#    scaleNote <- 'correlation scale'
#  }
  
  nc    <- 0
  vnam  <- matrix( unlist( strsplit(colnames(betaSig),'_') ),ncol=2,byrow=T)
  
  if( length(which(is.finite(match(colnames(y),vnam[,1])))) > 0 )nc <- 1
  if( length(which(is.finite(match(colnames(y),vnam[,2])))) > 0 )nc <- 2
  
  ix <- 1
  if(nc == 1)ix <- 2
  xnam <- vnam[,ix]
  vnam <- vnam[,nc]
  
  vnames <- unique(vnam)
  xnam <- unique(xnam[xnam != 'intercept'])
  
  if(SAVEPLOTS)pdf( file=.outFile(outfolder,'betaChains.pdf') ) # start plot
  
  mfrow <- .getPlotLayout(length(xnam))
  par(mfrow=mfrow, bty='n', oma=oma, mar=c(2,2,1,1), tcl= tcl, mgp=mgp)
  
  flist <- output$modelSummary$factorList
  if(length(flist) > 0){
    flist <- sort(unique(unlist(flist)))
  }
  
  for(k in 1:length(xnam)){
    
    tname <- xnam[k]
    
  #  kchain <- betaSig
  #  if(tname %in% flist){
  #    kchain <- bfSig
  #    tname <- rownames(eCont)[eCont[,tname] == 1]
  #  }
    
    tmp <- .chains2density(chains$bgibbs,varName=tname, cut=3)
    xt  <- tmp$x
    yt  <- tmp$y
    chainMat <- tmp$chainMat
    
    colF <- colorRampPalette(c('darkblue','orange'))
    cols <- colF(nrow(xt))
    
    snamek <- matrix( unlist(strsplit(colnames(chainMat),'_')),ncol=2,byrow=T)
    if(colnames(y)[[2]] %in% snamek[,1])nc <- 1
    if(colnames(y)[[2]] %in% snamek[,2])nc <- 2
    snamek <- snamek[,nc]
    
    nn <- nrow(chainMat)
    
    jk <- 1:ncol(chainMat)
    if(length(jk) > 20)jk <- sample(jk,20)
    plot(0,0,xlim=c(0,(1.4*nn)),ylim=range(chainMat[,jk]),
         xlab=' ',ylab=' ',cex=.01)
    
    for(j in jk){
      lines(chainMat[,j],col=cols[j])
      if(ncol(chainMat) < 15)text(nn,chainMat[nn,j],snamek[j],col=cols[j],pos=4)
      if(k == 1 & j == 1).plotLabel( paste('burn-in =',burnin),
                                     location='topright' )
    }
    .plotLabel(label=paste(letters[k],') ',tname,sep=''),location='topleft',above=T)
    
    abline(h=0,lwd=4,col='white')
    abline(h=0,lty=2)
    
    if(ncol(chainMat) >= 15) text(nn,mean(par('usr')[3:4]),
                                  paste(ncol(chainMat),'spp'),pos=4)
  }
  
  if(!SAVEPLOTS){
    readline('beta coefficient chains converged? -- return to continue ')
  } else {
    dev.off()
  }
  
  ######################### correlation chains, species at random
  
  
  if(SAVEPLOTS)pdf( file=.outFile(outfolder,'corChains.pdf') ) # start plot
  
    par(mfrow=c(2,2), bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
    
    w0 <- 1:ncol(chains$sgibbs)
    if(REDUCT){
      same <- .sameByColumn(chains$kgibbs)
      w0   <- order(same[lower.tri(same,diag=T)])
    } else {
      w0 <- sample(max(w0),80,replace=T)
    }
    ww <- 1:20
    
    for(jj in 1:4){
      
      ssj <- w0[ww]
      ssj <- ssj[is.finite(ssj)]
      if(length(ssj) == 0)break
      if(max(ssj) > max(w0))break
      ww  <- ww + 20
      
      tmp   <- .chains2density(rgibbsShort[,ssj])
      xt    <- tmp$x
      yt    <- tmp$y
      chainMat <- tmp$chainMat
      
      colF   <- colorRampPalette(c('black','brown','orange'))
      cols <- colF(nrow(xt))
      
      stk <- matrix( unlist(strsplit(colnames(chainMat),'_')),ncol=2,byrow=T)
      
      ws <- which(stk[,1] == stk[,2])
      if(length(ws) > 0){
        stk <- stk[-ws,]
        chainMat <- chainMat[,-ws]
      }
      
      rr <- range(chainMat)
      if(!is.finite(rr[1]) | !is.finite(rr[2]))next
      
      if(is.matrix(chainMat)){
        
        snamek <- stk[,1]
        nn <- nrow(chainMat)
        plot(0,0,xlim=c(0,(1.4*nn)),ylim=range(chainMat),xlab=' ',ylab=' ',cex=.01)
        
        jk <- 1:ncol(chainMat)
        if(length(jk) > 20)jk <- sample(jk,20)
        
        for(j in jk){
          lines(chainMat[,j],col=cols[j])
          if(ncol(chainMat) < 15)text(nn,chainMat[nn,j],snamek[j],col=cols[j],pos=4)
          if(j == 1).plotLabel( paste('burnin =',burnin),location='topright' )
        }
        
        abline(h=0,lwd=4,col='white')
        abline(h=0,lty=2)
        
        if(ncol(chainMat) >= 15) text(nn,mean(par('usr')[3:4]),
                                      paste(ncol(chainMat),'spp'),pos=4)
      }
    }
    
    if(!SAVEPLOTS){
      readline('correlation chains converged? -- return to continue ')
    } else {
      dev.off()
    }
  
  ############################### beta posteriors as boxes
    
    fnames <- names( output$parameterTables$fMu )
    
    if(length(bfSig) > 0){
      
      nc    <- 0
      vnam  <- matrix( unlist( strsplit(colnames(bfSig),'_') ),ncol=2,byrow=T)
      
      ix <- 1
      nc <- 2
      y1 <- which(vnam[,1] %in% snames)
      y2 <- which(vnam[,2] %in% snames)
      if(length(y1) > length(y2))ix <- 2
      if(ix == 2)nc <- 1
      
      # if( length(which(is.finite(match(colnames(y),vnam[,1])))) > 5 )nc <- 2
      # if( length( which(is.finite(match(colnames(y),vnam[,2]))) ) > 5 )nc <- 1
      
      # ix <- 1
      # if(nc == 1)ix <- 2
      xnam <- vnam[,ix]
      vnam <- vnam[,nc]
      
      xpNames <- .replaceString(fnames,':','X')
      xpNames <- .replaceString(xpNames,'I(','')
      xpNames <- .replaceString(xpNames,')','')
      xpNames <- .replaceString(xpNames,'^2','2')
      xpNames <- .replaceString(xpNames,'*','TIMES')
      
      for(j in 1:length(fnames)){
        
        wc <- which(xnam == fnames[j])
        if(length(wc) < 2)next
        
        plab <- paste('beta_',xpNames[j],'.pdf',sep='')
        if(SAVEPLOTS)pdf( file=.outFile(outfolder,plab) ) # start plot
        
        
        par(mfrow=c(1,1),bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
        
        if(length(wc) > 100)wc <- sample(wc,100)
        
        mat <- bfSig[,wc]
        ord <- order(colMeans(mat),decreasing=F)
        tnam <- vnam[ wc[ord] ]
        bx   <- boxCol[ match(tnam, snames) ]
        bb   <- specColor[ match(tnam, snames) ]
        ry   <- range(mat)
        ymin <- min(mat) - diff(ry)*.15
        ymax <- max(mat) + diff(ry)*.15
        
        tmp <- .boxplotQuant( mat[,ord],xaxt='n',outline=F,ylim=c(ymin,ymax),
                              col=bx, border=bb, xaxt='n',lty=1)
        abline(h=0,lwd=2,col='grey',lty=1)
        
        dy <- .05*diff(par()$yaxp[1:2])
        
        cext <- .fitText2Fig(tnam,fraction=1)
        text((1:length(wc)) - .1,dy + tmp$stats[5,],tnam,srt=70,pos=4,
             col=bb, cex=cext)
        
        pl    <- par('usr')
        xtext <- pl[1]
        ytext <- pl[3] + diff(pl[3:4])*.85
        #   text(xtext, ytext, expression(hat(bold(beta))),pos=4,cex=1.1)    
        #   .plotLabel(paste(fnames[j], scaleNote,sep=', '),location='topleft', 
        #              cex=1.0)   
        .plotLabel(fnames[j],location='topleft', cex=1.0)
        
        if(!SAVEPLOTS){
          readline('95% posterior -- return to continue ')
        } else {
          dev.off()
        }
      }
      
      #one figure
      
      if(length(fnames) > 1){
        
        if(SAVEPLOTS)pdf( file=.outFile(outfolder,'betaAll.pdf') )  
        
        npp <- length(which(table(match(xnam,fnames)) > 1))
        
        mfrow <- .getPlotLayout(npp)
        par( mfrow=mfrow, bty='n', oma=oma, mar=c(1,1,1,1), tcl= tcl, mgp=mgp )
        
        k <- 0
        for(j in 1:length(fnames)){
          
          wc <- which(xnam == fnames[j])
          if(length(wc) < 2)next
          
          k <- k + 1
          
          plab <- paste('beta_',xnames[j],'.pdf',sep='')
          
          if(length(wc) > 100)wc <- sample(wc,100)
          
          mat <- bfSig[,wc]
          ord <- order(apply(mat,2,mean),decreasing=F)
          tnam <- vnam[ wc[ord] ]
          bx   <- boxCol[ match(tnam, snames) ]
          bb   <- specColor[ match(tnam, snames) ]
          ry   <- range(mat)
          ymin <- min(mat) - diff(ry)*.15
          ymax <- max(mat) + diff(ry)*.15
          
          tmp <- .boxplotQuant( mat[,ord],xaxt='n',outline=F,ylim=c(ymin,ymax),
                                col=bx, border=bb, xaxt='n',lty=1)
          abline(h=0,lwd=2,col='grey',lty=1)
          
          dy <- .05*diff(par()$yaxp[1:2])
          
          cext <- .fitText2Fig(tnam,fraction=.9)
          
          text((1:length(wc)) - .1,dy + tmp$stats[5,],tnam,srt=70,pos=4,
               col=bb, cex=cext)
          
          pl    <- par('usr')
          xtext <- pl[1]
          ytext <- pl[3] + diff(pl[3:4])*.2
          #  text(xtext, ytext, expression(hat(bold(beta))),pos=4,cex=1.5)    
          .plotLabel(paste(letters[k],') ',fnames[j],sep=''),
                     location='topleft', cex=1.0)   
        }
        
        if(!SAVEPLOTS){
          readline('95% posterior -- return to continue ')
        } else {
          dev.off()
        }
      }
      
    }
    
  ############################### beta posteriors, traits
  
  if(TRAITS){
    
    M  <- nrow(specByTrait)
    nc    <- 0
    vnam  <- matrix( unlist( strsplit(colnames(chains$agibbs),'_') ),
                     ncol=2,byrow=T)
    mnames <- colnames(specByTrait)
    
    if( length(is.finite(match(mnames,vnam[,1]))) > 0 )nc <- 2
    if( length(is.finite(match(mnames,vnam[,2]))) > 0 )nc <- 1
    
    ix <- 1
    if(nc == 1)ix <- 2
    xnam <- vnam[,ix]
    vnam <- vnam[,nc]
    
    if(length(traitColor) == 1)traitColor <- rep(traitColor, M)
    tboxCol <- .getColor(traitColor,.4)
    
    traitSd <- apply(plotByTrait,2,sd,na.rm=T)
    traitSd <- matrix(traitSd,nrow(chains$agibbs),length(traitSd),byrow=T)
    
    for(j in 2:length(xnames)){
      
      wc <- which(xnam == xnames[j])
      if(length(wc) < 2)next
      
      if(SAVEPLOTS)pdf( file=.outFile(outfolder,'traits.pdf') ) # start plot
      
      par(mfrow=c(1,1),bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
      
      if(length(wc) > 100)wc <- sample(wc,100)
      
      mat <- chains$agibbs[,wc]*xSd[j]/traitSd
      ord <- order(apply(mat,2,mean),decreasing=F)
      tnam <- vnam[ wc[ord] ]
      bx   <- tboxCol[ match(tnam, mnames) ]
      bb   <- traitColor[ match(tnam, mnames) ]
      ry   <- range(mat)
      ymin <- min(mat) - diff(ry)*.15
      ymax <- max(mat) + diff(ry)*.15
      
      tmp <- .boxplotQuant( mat[,ord],xaxt='n',outline=F,ylim=c(ymin,ymax),
                            col='grey', border='black', xaxt='n',lty=1)
      abline(h=0,lwd=2,col='grey',lty=1)
      
      dy <- .05*diff(par()$yaxp[1:2])
      text((1:length(wc)) - .1,dy + tmp$stats[5,],tnam,srt=70,pos=4,
           col=bb, cex=.7)
      
      pl    <- par('usr')
      xtext <- pl[1]
      ytext <- pl[3] + diff(pl[3:4])*.2
   #   text(xtext, ytext, expression(hat(bold(alpha))),pos=4,cex=1.5)    
      .plotLabel(xnames[j],location='bottomleft')  
      
      if(!SAVEPLOTS){
        readline('95% posterior -- return to continue ')
      } else {
        dev.off()
      }
    }
  }
  
  
  ########### cluster analysis
  
  covx <- cov(x)
  covy <- cov(y[,notOmit])
  
  wo <- which(whichZero[,1] %in% other | whichZero[,2] %in% other)
  if(length(wo) > 0)whichZero <- whichZero[-wo,]
  wo <- which(whConZero[,1] %in% other | whConZero[,2] %in% other)
  if(length(wo) > 0)whConZero <- whConZero[-wo,]
  
  nsim <- 500
  if(S > 50)nsim  <- 100
  if(S > 100)nsim <- 20
  
  tmp <- eigen( ematrix[notOther,notOther] )
  
  eVecs   <- tmp$vectors
  eValues <- tmp$values
  rownames(eVecs) <- snames[notOther]
  
  
  
  
  if(!GRIDPLOTS){
    
    clusterIndex <- NULL
    clusterOrder <- NULL
    
    if(S >= 8){
      
      clusterDat <- .clusterPlot( dcor = ematrix , ncluster=ncluster, PLOT=F)
      colCode    <- clusterDat$colCode
      cord       <- rev(clusterDat$corder)
      dord       <- notOther[!notOther %in% omit][cord]
      
      clusterIndex <- clusterDat$clusterIndex
      clusterOrder <- clusterDat$corder
    }
 #   names(clusterOrder) <- names(clusterIndex)
    
 #   tmp <- .multivarEmat(bchains=chains$bgibbs[,keepBC],covx,snames,
 #                        orderB=notOther[cord], alpha=ematAlpha, nsim=nsim)
 #   Ematrix   <- tmp$bm
 #   whichZero <- tmp$whichZero
 #   whConZero <- tmp$whConZero
    
    return( list(summaryCoeffs = summaryCoeffs, fit = fit, ematrix = ematrix,
                 clusterIndex = clusterIndex, clusterOrder = clusterOrder) )
  }
  
  
  if(SAVEPLOTS)pdf( file=.outFile(outfolder,'clusterDataE.pdf') ) # start plot
  
  mag <- mar
  mag[4] <- max(mar[4],6)
  par(mfrow=c(1,2), cex=.7, oma=oma, mar=mag, tcl= tcl, mgp=mgp)
 
  LABELS <- T
  if(S > 100 | !specLabs)LABELS <- F
  
  dcor <- .cov2Cor(covy)
  dcor[is.na(dcor)] <- 0
  
  tmp <- .clusterPlot( dcor =  dcor,main='',cex=.2,
                       ncluster=ncluster,colCode=specColor[notOmit], 
                       textSize=.4, LABELS = LABELS)
  colCode <- tmp$colCode
  
  clusterIndex <- tmp$clusterIndex
  clusterOrder <- tmp$corder
  
  .plotLabel('a) Data correlation',above=T, cex=1.7)
  
  tmp <- .clusterPlot( dcor = ematrix[notOther,notOther] ,main='',cex=.2,
                       ncluster=ncluster, 
                       colCode=colCode, textSize=.5, LABELS = LABELS)
  .plotLabel('b) E correlation',above=T, cex=1.7)
  
  clusterIndex <- cbind( clusterIndex, tmp$clusterIndex )
  clusterOrder <- cbind( clusterOrder, tmp$corder )
  
  rownames(clusterIndex) <- rownames(clusterOrder) <- snames[notOmit]
  colnames(clusterIndex) <- colnames(clusterOrder) <- c('data','E')
  
  if(!SAVEPLOTS){
    readline('Data and E responses to X -- return to continue ')
  } else {
    dev.off()
  }
  
  ########### E communities
  
  imat <- output$y
  imat[imat > 0] <- 1
  iord <- colSums(imat)
  etab  <- table(clusterIndex[,'E'])
  eComs <- matrix(NA,ncluster, max(etab))
  ename <- rep( character(0), max(etab) )
  
  egroup <- clusterIndex[,'E']
  bTab   <- cbind(egroup,bTab[notOther,])
  summaryCoeffs$betaSig <- bTab
  
  bfTab <- cbind(egroup, bfTab[notOther,])
  summaryCoeffs$fBetaSig <- bfTab
  
  for(j in 1:ncluster){
    
    wj <- which(clusterIndex[,'E'] == j)
    jname <- rownames(clusterIndex)[wj]
    jname <- jname[order(iord[jname],decreasing=T)]
    eComs[j,1:length(jname)] <- jname
    mm    <- min( c(3,length(jname)) )
    jj    <- substr(jname[1:mm],1,6)
    ename[j] <- paste0(jj,collapse='_')
  } 
  rownames(eComs) <- ename
  eComs <- t(eComs)
  
  ########### ordination
  
  if(SAVEPLOTS)pdf( file=.outFile(outfolder,'ordination.pdf') ) # start plot
  
  clusNames <- eComs[1,]
  
  lambda <- eValues/sum(eValues)
  cl     <- cumsum(lambda)
  
  cbord <- .getColor(specColor[notOther],.4)
  cfill <- .getColor(specColor[notOther],.4)
  
  par(mfcol=c(2,2), bty='n', cex = cex, mar=c(4,4,1,1))
  
  p1 <- paste('Axis I (',round(100*lambda[1],0),'%)',sep='')
  p2 <- paste('Axis II (',round(100*lambda[2],0),'%)',sep='')
  p3 <- paste('Axis III (',round(100*lambda[3],0),'%)',sep='')
  
  xlim <- range(eVecs[,1])
  
  plot(eVecs[,1],eVecs[,2],cex=1,col=cbord, bg = cfill, pch=16,
       xlab=p1, ylab = p2) 
  abline(h=0,col=.getColor('black',.1),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.1),lwd=2,lty=2)
  
  text(eVecs[clusNames,1],eVecs[clusNames,2],substr(clusNames,1,7))
  
  plot(eVecs[,1],eVecs[,3],cex=1,col=cbord, bg = cfill, pch=16,
       xlab=p1, ylab = p3) 
  abline(h=0,col=.getColor('black',.1),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.1),lwd=2,lty=2)
  
  text(eVecs[clusNames,1],eVecs[clusNames,3],substr(clusNames,1,7))
  
  plot(eVecs[,2],eVecs[,3],cex=1,col=cbord, bg = cfill, pch=16,
       xlab=p2, ylab = p3)
  abline(h=0,col=.getColor('black',.1),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.1),lwd=2,lty=2)
  
  text(eVecs[clusNames,2],eVecs[clusNames,3],substr(clusNames,1,7))
  
  plot(cl,type='s',xlab='Rank',ylab='Proportion of variance',xlim=c(.9,S),
       ylim=c(0,1),log='x')
  lines(c(.9,1),c(0,cl[1]),lwd=2,type='s')
  for(j in 1:length(lambda))lines(c(j,j),c(0,cl[j]),col='grey')
  lines(cl,lwd=2,type='s')
  abline(h=1,lwd=2,col=.getColor('grey',.5),lty=2)
  
  if(!SAVEPLOTS){
    readline('ordination of E matrix -- return to continue ')
  } else {
    dev.off()
  }
  
  ########### dimension reduction ############
  
  if(REDUCT){ 
    
    graphics.off()
    
    if(SAVEPLOTS)pdf( file=.outFile(outfolder,'dimRed.pdf') ) # start plot
    
    mk <- .modalValuesInArray(chains$kgibbs,2)[notOmit]
    NK <- table( table(mk) )
    mk <- length(NK)
    
    r <- otherpar$r
    
    par(bty='n')
    scale <- SO/3
    if(SMALLPLOTS)scale <- 10*scale
    .mapSetup(c(1,SO),c(1,SO),scale=scale)
    
    xl <- SO/15
    yl <- SO/8
    
    en <- SO*(SO+1)/2
    
    plot(0,0,xlim=c(0,SO+xl),ylim=c(0,SO+xl),cex=.01,xaxt='n',yaxt='n',
         xlab=' ',ylab=' ')
    
    rect(xl,yl,SO+xl,SO+yl,col='wheat',border='wheat',lty=2,lwd=2)
    polygon(c(xl,SO+xl,xl),c(yl,yl,SO+yl),col='blue',border='darkblue')
    rect(0,yl/10,r,mk+yl/10,col='turquoise',border='brown')
    arrows(r, yl/20 + mk/3, r + SO/5, yl/20*(mk + 1), col='brown',
           angle=20, lwd=2)
    text(xl+SO/4,yl+SO/3,bquote(Sigma == .(en)), col='wheat', cex=1.4 )
    text(r + SO/5, yl/20*(mk + 1),
         paste('Z (',mk,' x ',r,' = ',mk*r,')',sep=''),col='brown',
         cex=1.,pos=4)
    .plotLabel('Dimensions',above=T)
    
    if(!SAVEPLOTS){
      readline('reduction from sigma to Z -- return to continue ')
    } else {
      dev.off()
    }
  } 
  
  ########### grid/correlation analysis
  
  
  if(SAVEPLOTS)pdf( file=.outFile(outfolder,'clusterGridR.pdf') )
  
  spl <- specLabs
  if(S > 80)spl <- F
  
  par(mfrow=c(1,1),bty='n',cex=1, oma=oma, mar=mag, tcl= tcl, mgp=mgp)
  
  colnames(corMu) <- rownames(corMu) <- colnames(y)
  
  psize <- .62
  if(SMALLPLOTS)psize <- psize/2
  
  par(plt=c(.03,.15,.1,.9), bty='n', new=F)
  tmp <- .clusterPlot( dcor = corMu[notOmit,notOmit] ,main=' ',cex=.2,
                       ncluster=ncluster, 
                       colCode=specColor[notOmit], textSize=.5, 
                       LABELS = F)
  colCode   <- tmp$colCode
  corder    <- rev(tmp$corder)
  specOrder <- snames[notOmit[corder]]
  
  clusterIndex <- cbind( clusterIndex, tmp$clusterIndex )
  clusterOrder <- cbind( clusterOrder, tmp$corder )
  
  ncc <- ncol(clusterIndex)
  colnames(clusterIndex)[ncc] <- colnames(clusterOrder)[ncc] <- 'R'
  
  if(LABELS){
    
    par(plt=c(.15,.33,.1,.9), bty='n', new=T)
    plot(c(0,0),c(0,0),col='white',xlim=range(c(0,1)),ylim=c(0,SO),
         xaxt='n',yaxt='n',xlab='',ylab='')
    xl <- rep(.5,SO)
    
    yl <- c(1:SO) + par('usr')[3] - .5
    cex <- .fitText2Fig(specOrder,fraction=1.2)
    text( xl,yl,rev(specOrder),pos=3,cex=cex, col=rev(colCode[corder]))
  }
  
  knames <- snames[notOmit]
  
  sgibbs <- chains$sgibbs
  if(REDUCT)sgibbs <- output$chains$sgibbs
  
  tmp <- .invMatZero(sgibbs,nsim=nrow(sgibbs),snames=snames,
                     knames=specOrder,index=NULL, COMPRESS=T, 
                     REDUCT=REDUCT,
                     sigErrGibbs = output$chains$sigErrGibbs, 
                     kgibbs = output$chains$kgibbs,
                     otherpar = otherpar, alpha=ematAlpha)
  marIn <- tmp$inMarMat
  conIn <- tmp$inConMat
  wm    <- which(marIn[,1] %in% omit |  marIn[,2] %in% omit)
  if(length(wm) > 0)marIn <- marIn[-wm,]
  wm    <- which(conIn[,1] %in% omit |  conIn[,2] %in% omit)
  if(length(wm) > 0)conIn <- conIn[-wm,]
  
  sigCor <- c(nrow(marIn),nrow(conIn))/SM/(SM - 1)
  sigCor <- round(100*sigCor,0)
  names(sigCor) <- c('n_marIn','n_conIn')
  
  mor <- notOmit[corder]
  
  crr <- corMu[mor,mor]
  marIn[,1] <- match(marIn[,1],mor)
  marIn[,2] <- match(marIn[,2],mor)
  conIn[,1] <- match(conIn[,1],mor)
  conIn[,2] <- match(conIn[,2],mor)
  
  makeCR <- list('white' = conIn,'grey' = marIn)
  if(!is.null(specColor))textCol = colCode[mor]
  
  par(plt=c(.33, .33 + psize,.1,.9), bty='n', new=T)
  
  slim <- quantile(crr[lower.tri(crr)],c(.05,.95))
  
  .corPlot(crr, slim=slim, makeColor=makeCR,plotScale=.99,
           PDIAG=T,corLines=corLines, textCol = colCode[corder], 
           specLabs = F, squarePlot = F,
           textSize=1, widex = width, widey = height, new=T, add=F)
  
  ll <- paste(c('Cond Ind (white) = ', 'Cond & Marg Ind (grey) = '),
              sigCor,c('%','%'),sep='')
  legend('topright',ll,bty='n',cex=.8)
  
  .plotLabel(expression( paste(hat(bold(R)),'structure'  )),above=T, cex=.9)
  
  if(!SAVEPLOTS){
    readline('posterior correlation for model -- return to continue ')
  } else {
    dev.off()
  }
  
  ########################### cluster Fmat with beta
  
  if(Q > 4){
    
    graphics.off()
    if(SAVEPLOTS)pdf( file=.outFile(outfolder,'gridF_B.pdf') ) # start plot
    
    main1 <- expression( paste('Sensitivity ',hat(F)))
    main2 <- expression( paste('Responses ',hat(B)))
    
    fBetaMu <- output$parameterTables$fBetaMu
    
    mat1 <- .cov2Cor(fMat)
    mat2 <- fBetaMu
    .clusterWithGrid(mat1, mat2, expand=5, mainLeft=main1, DIST=F,
                     main1=main1, main2 = main2,
                     leftClus=T, rightClus=F, topClus2=T, 
                     leftLab=F, rightLab=T, topLab1=T, topLab2=F,
                     colOrder1=NULL, colOrder2=NULL, rowOrder=NULL, 
                     colCode1 = NULL, colCode2 = boxCol[notOther],
                     rowCode=NULL,lower1 = T, diag1 = T,lower2 = F, 
                     diag2 = F, slim1=NULL, slim2=NULL)
    
    if(!SAVEPLOTS){
      readline('F & beta structure -- return to continue ')
    } else {
      dev.off()
    } 
  }
  #################################### cluster Emat
  
  graphics.off()
  if(SAVEPLOTS)pdf( file=.outFile(outfolder,'clusterGridE.pdf') ) # start plot
  
  mat1 <- ematrix[notOther,notOther]
  main1 <- expression(paste('Species ',hat(E)))
  .clusterWithGrid(mat1, mat2=NULL, expand=1, mainLeft=main1, DIST=F,
                   main1=NULL, main2 = NULL,
                   leftClus=T, rightClus=F, topClus1=F, topClus2=F, 
                   leftLab=T, rightLab=F, topLab2=F,
                   colOrder1=NULL, colOrder2=NULL, rowOrder=NULL, 
                   colCode1 = boxCol[notOther],
                   rowCode=boxCol[notOther],
                   lower1 = T, diag1 = T,lower2 = F, diag2 = F,
                   slim1=NULL, slim2=NULL, horiz1=clusterIndex[,'E'])
  
  if(!SAVEPLOTS){
    readline('E: model-based response to X -- return to continue ')
  } else {
    dev.off()
  }
  
  ################# resid and Egrid
  
  graphics.off()
  if(SAVEPLOTS)pdf( file=.outFile(outfolder,'gridR_E.pdf') ) # start plot
  
  mat1 <- crr
  mat2 <- ematrix[notOther,notOther]
  main1 <- expression(paste('Ordered by error ',hat(R)))
  main2 <- expression(paste('Response ',hat(E)))
  
  .clusterWithGrid(mat1, mat2=mat2, expand=1, mainLeft='Species', DIST=F,
                   main1=main1, main2 = main2,
                   leftClus=T, rightClus=F, topClus1=F, topClus2=F, 
                   leftLab=T, rightLab=F, topLab2=F,
                   colOrder1=NULL, colOrder2=NULL, 
                   rowOrder=NULL, 
                   colCode1 = NULL,
                   rowCode = NULL,
                   lower1 = T, diag1 = T,lower2 = T, diag2 = T,
                   slim1=NULL, slim2=NULL)
  # rowCode=rev(colCode[corder][notOther]),
  if(!SAVEPLOTS){
    readline('comparison R vs E -- return to continue ')
  } else {
    dev.off()
  }
  
  ################# data vs E grid
  
  graphics.off()
  if(SAVEPLOTS)pdf( file=.outFile(outfolder,'gridY_E.pdf') ) # start plot
  
  ytmp <- jitter(y[,mor],1e-10)
  cory <- cor(ytmp)
  
  mat1 <- cory
  mat2 <- ematrix[notOther,notOther]
  main1 <- 'Ordered by data, cor(Y)'
  main2 <- expression(paste('Response ',hat(E)))
  
  .clusterWithGrid(mat1, mat2=mat2, expand=1, mainLeft='Species', DIST=F,
                   main1=main1, main2 = main2,
                   leftClus=T, rightClus=F, topClus1=F, topClus2=F, 
                   leftLab=T, rightLab=F, topLab2=F,
                   colOrder1=NULL, colOrder2=NULL, rowOrder=NULL, 
                   colCode1 = NULL,
                   rowCode=NULL,
                   lower1 = T, diag1 = T,lower2 = T, diag2 = T,
                   slim1=NULL, slim2=NULL)
  
  if(!SAVEPLOTS){
    readline('raw data vs E -- return to continue ')
  } else {
    dev.off()
  }
  
  #################### beta grid
  if(betaGrid & nrow(fBetaMu) > 4){
    
    graphics.off()
    
    if(SAVEPLOTS)pdf( file=.outFile(outfolder,'clusterGridB.pdf') ) # start plot
    
    mat1 <- ematrix[notOther,notOther]
    mat2 <- t(fBetaMu)
    main1 <- expression(paste('Species ',hat(E)))
    main2 <- expression(paste(hat(B),' by predictor'))
    
    ee <- ncol(mat2)/(ncol(mat1) + ncol(mat2) )
    ee <- max(ee,.05)
    .clusterWithGrid(mat1, mat2, expand=ee, mainLeft=main1, DIST=F,
                     main1=main1, main2 = main2,
                     leftClus=F, rightClus=F, topClus1=T, topClus2=T, 
                     leftLab=F, rightLab=F, topLab2=T,
                     colOrder1=NULL, colOrder2=NULL, rowOrder=NULL, 
                     colCode1 = boxCol[notOther],
                     rowCode=boxCol[notOther],
                     lower1 = T, diag1 = T,lower2 = F, diag2 = F,
                     slim1=NULL, slim2=NULL, vert1=clusterIndex[,'E'],
                     horiz2=clusterIndex[,'E'])
    
    if(!SAVEPLOTS){
      readline('beta ordered by response to X -- return to continue ')
    } else {
      dev.off()
    }
  
     if(TRAITS){

       bb <- betaTraitMu[-1,]
       ord <- order(colSums(abs(bb)),decreasing=T)
       bb  <- bb[,ord]
       bl  <- bb[,ord]
       bh  <- bb[,ord]
       ror <- order(rowSums(abs(bb)),decreasing=T)
       bb  <- bb[ror,]
       bl  <- bl[ror,]
       bh  <- bh[ror,]
       
       white <- which(bl < 0 & bh > 0,arr.ind=T)
       
       makeColor <- list('white' = white )
       
       if(SAVEPLOTS)pdf( file=.outFile(outfolder,'gridTraitB.pdf') ) 
       
       plotScale <- max(c(10,c(S,Q)/10))
       
       par(mfrow=c(1,1), bty='n', oma=c(1,1,1,1), 
           mar=c(5,4,4,2), tcl= tcl, mgp=mgp)
       
       ht <- nrow(bb)/ncol(bb)*width
       
       .corPlot(bb,slim=NULL,PDIAG=F,plotScale=plotScale,
                makeColor=NULL,textSize=.7, specLabs = specLabs,
                corLines=corLines,tri='both', 
                cex=.6,squarePlot=F,LEGEND=T,
                widex=width,widey=ht,new=F)
       
       if(!SAVEPLOTS){
         readline('trait beta -- return to continue ')
       } else {
         dev.off()
       }
     }
  }
  
  all <- list(summaryCoeffs = summaryCoeffs, fit = fit, ematrix = ematrix,
              eComs = eComs,
              clusterIndex = clusterIndex, clusterOrder = clusterOrder,
              eVecs = eVecs, eValues = eValues) 
  all <- all[ order(names(all)) ]
}
  
  
.gjamPrediction <- function(output, newdata, y2plot, PLOT, ylim){
  
  xdata <- ydataCond <- NULL
  tiny  <- 1e-10
  wHold <- phiHold <- ploHold <- sampleWhold <- NULL
  
  ng     <- output$modelList$ng
  burnin <- output$modelList$burnin
  
  nsim <- 500
  if('nsim' %in% names(newdata))nsim <- newdata$nsim
  
  if( is.null(newdata) ){
    
    if(PLOT){
      
      y1 <- output$y
      y2 <- output$modelSummary$ypredMu
      if(!is.null(y2plot)){
        y1 <- y1[,y2plot]
        y2 <- y2[,y2plot]
      }
      
      tmp <- .bins4data(y1)
      breaks <- tmp$breaks
      bins   <- tmp$bins
      nbin   <- tmp$nbin
      
      if(length(bins) > 0){
        breaks <- bins
        nPerBin <- NULL
      }
      
      .plotObsPred(y1,y2, nPerBin = NULL, breaks=breaks, ylimit = ylim,
                   fill='lightblue', box.col='darkblue',POINTS=F)
      abline(0,1,lwd=4,col='white')
      abline(0,1,lwd=2,col='grey',lty=2)
    }
    return(  list( ypredMu = output$modelSummary$ypredMu, 
                   ypredSe = output$modelSummary$ypredSd ) )
  }
  
  S <- SO <- S1 <- ncol(output$y)
  Q <- ncol(output$x)
  n <- nrow(output$x)
  y <- output$y
  x <- output$x
  
  xnames <- colnames(x)
  ynames <- colnames(y)
  
  cindex <- NULL
  
  notOther <- c(1:S)
  other    <- grep('other',ynames)
  if(length(other) > 0)notOther <- notOther[-other]
  SO       <- length(notOther)
  
  otherpar <- output$otherpar
  censor   <- output$censor
  REDUCT   <- output$REDUCT
  
  notStandard <- output$modelList$notStandard
  
  
  NEWX <- F
  if('xdata' %in% names(newdata))NEWX <- T
  
  
  effort <- output$effort
  effMat <- output$y*0 + 1
  if(!is.null(effort)){
  }
  
  if('effort' %in% names(newdata)){
    effort <- newdata$effort
    effMat <- matrix(1,nrow(newdata$xdata),ncol(output$y))
    effMat[,effort$columns] <- effort$values
  }
  effort <- list(columns = c(1:S), values = effMat)

  if(REDUCT){
    N  <- output$otherpar$N
    r  <- output$otherpar$r
  }
  cuts <- output$parameterTables$cutMu
  cuts <- cbind(-Inf,0,cuts,Inf)
  
  isFactor   <- output$modelSummary$isFactor
  factorList <- output$modelSummary$factorList
  contrasts  <- output$modelSummary$contrasts
  formula    <- output$modelList$formula
  xscale     <- output$standX
  if(is.matrix(xscale))  xscale <- t(xscale)
  facNames   <- names(factorList)
  
  typeNames <- output$modelSummary$typeNames
  tmp       <- .gjamGetTypes(output$modelSummary$typeNames)
  typeFull  <- tmp$typeFull
  typeCols  <- tmp$typeCols
  allTypes  <- unique(typeCols)
  typeCode  <- tmp$TYPES[typeCols]
  FCgroups  <- attr(typeNames,'FCgroups')
  CCgroups  <- attr(typeNames,'CCgroups')
  CATgroups <- attr(typeNames,'CATgroups')
  
  xdataNames <- output$modelSummary$xdataNames
  ydataNames <- output$modelSummary$ydataNames
  
  COND     <- T
  condCols <- numeric(0)
  
  if( NEWX ){   ################ out-of-sample
    
    xdata <- newdata$xdata
    COND  <- F
    nx    <- n <- nrow(newdata$xdata)
    
    effMat <- matrix(1, nx, S)
    
    if( 'effort' %in% names(newdata) ){
      ev     <- newdata$effort$values
      effMat <-  matrix(1, nx, S)
      effMat[,newdata$effort$columns] <- ev
    }
    effort <- list(columns = c(1:S), values = effMat)

    if( 'ydataCond' %in% names(newdata) )
      stop('supply either xdata or ydataCond, not both')
    
    ydataCond <- NULL
 #   n <- nrow(xdata)
    
    if(length(factorList) > 0){
      
      for(j in 1:length(factorList)){
        
        nf <- names(factorList)[j]
        wf <- which(names(xdata) == nf)
        wo <- which(names(output$xdata) == nf)
        wc <- which(names(output$modelSummary$contrasts) == names(factorList)[j])
        cc <- output$modelSummary$contrasts[[wc]]
        
        xdata[[wf]] <- factor( xdata[[wf]], levels = levels(output$xdata[[wo]]) )
        
        attr(xdata[[wf]],'contrasts') <- cc
     #   attr(xdata[[wf]],'reference') <- rownames(cc)[!rownames(cc) %in% colnames(cc)]
      }
    }
  
    x <- matrix(0,nx,Q)
    y <- matrix(0,nx,S)
    colnames(x) <- xnames
    colnames(y) <- ynames
    yp <- y
    
    standRows <- output$modelSummar$standRows
    standMat  <- output$modelSummar$standMat
    
    tmp <- .gjamXY(formula, xdata, yp, typeNames, 
                   notStandard, checkX = F, xscale = xscale)
    
    x  <- tmp$x
    ws <- which(xscale[1,] != 0)
    if(length(ws) > 0){
      xx <- xdata[,colnames(xscale)[ws],drop=F]
      for(k in 1:length(ws)){
        x[,ws[k]] <- (xx[,k] - xscale[1,ws[k]])/xscale[2,ws[k]]
      }
    }

    yp <- w <- x%*%output$parameterTables$betaMu
    
    wca <- which(typeNames == 'CA')
    if(length(wca) > 0){
      yp[,wca][yp[,wca] < 0] <- 0
    }
    
    wda <- which(typeNames == 'DA')
    if(length(wda) > 0){
      yp[,wda] <- round(yp[,wda]*effMat[,wda],0)
      yp[,wda][yp[,wda] < 0] <- 0
    }

    ordCols <- which(typeNames == 'OC')
    if(length(ordCols) > 0){
      
      tmp   <- .gjamGetCuts(yp + 1,ordCols)
      cutLo <- tmp$cutLo
      cutHi <- tmp$cutHi
      
      for(k in ordCols){
        yp[,k] <- findInterval(yp[,k],cuts[k,]) - 1
      }
    }
    
    if(length(FCgroups) > 0){
      ntt <- max(FCgroups)
      for(i in 1:ntt){   
        wk      <- which( FCgroups == i )
        wo      <- which(wk %in% notOther)
        yp[,wk] <- .gjamCompW2Y(yp[,wk], notOther=wo)$ww
      }
    }
    
    if(length(CCgroups) > 0){
      ysum <- rep(1000,n)                   # CC use sum of 1000
      ntt  <- max(CCgroups)
      if(ntt > 0){
      for(i in 1:ntt){  ## normalize y 
        wk      <- which( CCgroups == i )
        wo      <- which(wk %in% notOther)
        yp[,wk] <- .gjamCompW2Y(yp[,wk], notOther=wo)$ww
        yp[,wk][yp[,wk] < 0] <- 0
        yp[,wk] <- round( sweep(yp[,wk],1,ysum,'*'), 0) 
      }
      }
    }
    
    tmp <- .gjamSetup(typeNames, x, yp, breakList=NULL, holdoutN=NULL, 
                      holdoutIndex=NULL,censor=NULL, effort=effort) 
 #   w <- tmp$w; z <- tmp$z; yp <- tmp$y; 
    other <- tmp$other
 #   plo      <- tmp$plo; phi <- tmp$phi
    ordCols  <- tmp$ordCols; disCols <- tmp$disCols; compCols <- tmp$compCols 
    minOrd   <- tmp$minOrd;   maxOrd <- tmp$maxOrd;  censorCA <- tmp$censorCA
    censorDA <- tmp$censorDA;   ncut <- ncol(cuts);   corCols <- tmp$corCols
    catCols  <- which(attr(typeNames,'CATgroups') > 0)
    sampleW  <- tmp$sampleW*0 + 1
    
    byCol <- byRow <- F
    if(attr(sampleW,'type') == 'cols')byCol <- T
    if(attr(sampleW,'type') == 'rows')byRow <- T
    indexW <- attr(sampleW,'index')
    
    
 #   tmp <- .gjamSetup(typeNames, x, yp*0, breakList=NULL, holdoutN=NULL, 
 #                     holdoutIndex=NULL,censor=NULL, effort=effort) 
 #   plo      <- tmp$plo
    
    pmax <- apply(output$y/output$effort$values,2,max) #max w at fitted effort
    
    phi  <- 2*matrix(pmax,n,S,byrow=T)                 #max w
    
    phi[,ordCols]  <- length(ordCols) + 10
    phi[,compCols] <- 2
    phi[,catCols]  <- 10
    
    plo <- -phi
    
    byCol <- byRow <- F
    if(attr(sampleW,'type') == 'cols')byCol <- T
    if(attr(sampleW,'type') == 'rows')byRow <- T
    indexW <- attr(sampleW,'index')
    
    cdex <- c(1:S)

  } else {                          ############## conditional
    
    if( !'ydataCond' %in% names(newdata) )
      stop('supply either xdata or ydataCond, not both')
    xdata     <- output$xdata
    ydataCond <- newdata$ydataCond
    condNames <- colnames(ydataCond)
    
    if('other' %in% condNames){
      condNames <- condNames[condNames != 'other']
      ydataCond <- ydataCond[,condNames]
    }
    
    x <- output$x
    y <- yp <- output$y 
    n <- nrow(x)
    
    condCols <- which(colnames(yp) %in% condNames)
    
    colnames(ydataCond) <- .replaceString(colnames(ydataCond),'-','')
    colnames(ydataCond) <- .replaceString(colnames(ydataCond),'_','')
    colnames(ydataCond) <- .replaceString(colnames(ydataCond),' ','')
    
    yp[,colnames(ydataCond)] <- ydataCond
    
    tmp <- .gjamSetup(typeNames, x, yp, breakList=NULL, holdoutN=NULL, 
                      holdoutIndex=NULL,censor=NULL, effort=effort) 
    w <- tmp$w; z <- tmp$z; yp <- tmp$y; other <- tmp$other
    plo      <- tmp$plo; phi <- tmp$phi
    ordCols  <- tmp$ordCols; disCols <- tmp$disCols; compCols <- tmp$compCols 
    minOrd   <- tmp$minOrd;   maxOrd <- tmp$maxOrd;  censorCA <- tmp$censorCA
    censorDA <- tmp$censorDA;   ncut <- ncol(cuts);   corCols <- tmp$corCols
    effort   <- tmp$effort
    catCols  <- which(attr(typeNames,'CATgroups') > 0)
    sampleW  <- tmp$sampleW
    sampleW[,-condCols] <- 1
    
    
    standMat <-  output$modelSummary$standMat
    standRows <- output$modelSummary$standRows
    
    tmp <- .gjamXY(formula, xdata, yp, typeNames, 
                   output$modelList$notStandard, checkX = F, 
                   xscale = xscale)
    xdata <- tmp$xdata
    standRows <- tmp$standRows
    standMat  <- tmp$standMat;    standMu <- tmp$standMu
    
    byCol <- byRow <- F
    if(attr(sampleW,'type') == 'cols')byCol <- T
    if(attr(sampleW,'type') == 'rows')byRow <- T
    indexW <- attr(sampleW,'index')
    
    wy   <- match( colnames(ydataCond), colnames(y) )
    ny   <- max(wy)    # or length?!
    if(is.na(ny))stop( 'colnames(newdata$ydataCond) must occur in output$y' )
    
    cdex <- c(1:S)[-wy]
    
    CCsums <- numeric(0)
    if(!is.null(CCgroups)){
      ncc    <- max(CCgroups)
      for(j in 1:ncc){
        wjk    <- which(CCgroups == j)
        CCsums <- append(CCsums,list( rowSums(y[,wjk]) ) )
      }
    }
  } ##############################
  
  
  if(length(other) > 0)cdex <- cdex[!cdex %in% other]
  S1   <- length(cdex)
  
 # yg <- yp <- y
  
  yg <- yp
  
  FULL <- F
  if(length(yp) < 10000) FULL <- T
  
  if(FULL){
    ygibbs <- wgibbs <- matrix(0,nsim,length(yp))
  }
  
 # wmax <- ymax <- apply(y,2,max)
 # wmax <- wmax/effMat
   
  .updateW <- .wWrapper(REDUCT, n, S, effMat, corCols, typeNames, 
                        typeFull, typeCols, 
                        allTypes, holdoutN=0, holdoutIndex=NULL, censor, 
                        censorCA, censorDA, notOther, sampleW, byRow, byCol,
                        indexW, ploHold, phiHold, sampleWhold)
  sigmaerror <- otherpar$sigmaerror
  
  
  ypred  <-  ypred2 <- wcred <- wcred2 <- matrix(0,n,S)
  gvals  <- sample(burnin:ng,nsim,replace=T)
  
  pbar <- txtProgressBar(min=1,max=nsim,style=1)
  ig   <- 0
  
  corColC <- cdex[cdex %in% corCols]
  corColW <- which(cdex %in% corCols)
  
  ddex <- which(notOther %in% cdex)
  
  cutg <- cuts
  ncut <- ncol(cutg)
  
  ccols <- which(typeNames != 'CON')
  
  kg     <- 1
  rndEff <- w*0
  
  prPresent <- w*0
  
  ############ E matrix
  emat <- matrix(0,S,S)
  colnames(emat) <- rownames(emat) <- ynames
  lo <- hi <- lm <- hm <- ess <- emat
  
  eCont <- output$modelSummary$eCont
  dCont <- output$modelSummary$dCont
  
  covE <- cov( x%*%dCont )  # note that x is standardized
  
  
  nfact <- length(output$modelSummary$factorList)
  eCont <- output$modelSummary$eCont
  frow  <- NULL
  if(nfact > 0){
    frow <- rep(0,Q)
    for(j in 1:nfact){
      frow[ match(factorList[[j]], xnames) ] <- j
    }
  }
      
  
  ############
  
  wxx <- which(notStandard %in% xnames)
  nss <- NULL
  if(length(wxx) > 0){
    nss <- match(notStandard[wxx],xnames)
  }

  lCont <- output$modelSummary$eCont
  
  q1 <- nrow(eCont)
  fnames   <- rownames(eCont)
  facList2 <- factorList
  if(nfact > 0){
    for(j in 1:nfact){
      wj <- which(names(xdata) == names(factorList)[j])
      facList2[[j]] <- levels(xdata[[wj]])
    }
  }
  
  
  
  for(g in gvals){
    
    bg   <- matrix( output$chains$bgibbs[g,], Q, S)
    bg[standRows,] <- bg[standRows,]*standMat[standRows,]
    
    muw  <- x%*%bg
    
    if(REDUCT){
      Z  <- matrix(output$chains$sgibbs[g,],N,r)
      sg <- .expandSigma(output$chains$sigErrGibbs[g], S, Z = Z, 
                         output$chains$kgibbs[g,], REDUCT = T)
      K    <- output$chains$kgibbs[g,]
      covR <- .solveArma( (1/sigmaerror)*crossprod(Z[K,]) + diag(r) ) # Sigma_W
      z1   <- crossprod( Z[K,]/sigmaerror,t(w - muw) )        
      RR   <- .rmvnormArma(n, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
      endEff <- RR%*%t(Z[K,])
    } else {
      sg <- .expandSigma(output$chains$sgibbs[g,], S = S, REDUCT = F)
    }
    
    alpha <- .sqrtRootMatrix(bg,sg,DIVIDE=T)
    
    #################
    
    bgg <- bg[,notOther]
    if(!is.null(nss))bgg[nss,] <- bgg[nss,]*standMat[nss,notOther]
    
    agg <- .sqrtRootMatrix(bgg,sg[notOther,notOther],DIVIDE=T)  #cor-stand scale
    
    if(nfact > 0){
      agg <- lCont%*%agg    #standardized for x and cor scale for y
      for(k in 1:nfact){
        fk  <- facList2[[k]]
        mua <- colMeans(agg[drop=F,fk,])
        nl  <- length(fk)
        agg[fk,] <- agg[fk,] - matrix(mua,nl,SO,byrow=T)
      }
    } else {
      agg <- agg[drop=F,-1,]
    }
    
    egg         <- lCont%*%bgg          #standardized for x, not cor for y
 #   fsens       <- egg%*%sinv%*%t(egg)
    ###################
    
  
    if( 'OC' %in% typeCode ){
      cutg[,3:(ncut-1)] <- matrix( output$chains$cgibbs[g,], S)
      tmp   <- .gjamGetCuts(yg + 1,ordCols)
      cutLo <- tmp$cutLo
      cutHi <- tmp$cutHi
      
      plo[,ordCols] <- cutg[cutLo]
      phi[,ordCols] <- cutg[cutHi]
    }
    
    tmp <- .updateW( x, w, yg, muw, sg, alpha, cutg, plo, phi, 
                     rndEff=rndEff, sigmaerror, wHold )
    w   <- tmp$w
    if( !COND )yg  <- tmp$yp   
      
  #  plo <- tmp$plo   # changes only for CAT
  #  phi <- tmp$phi
    
    
    if(COND){ 
      tmp <- .conditionalMVNRcpp(w, muw, sg, cdex = ddex, S)  
      muc    <- tmp$mu
      sgp    <- tmp$vr
      if(S1 == 1){
        w[,ddex] <- matrix(rnorm(n,muc,sqrt(sgp[1])))
      } else {
        w[,ddex] <- .rMVN(n,muc,sgp)
      }
      muw[,ddex] <- muc
  #    yg[,ddex] <- w[,ddex]
      
      if( length(corColC) > 0 ){    #expanded w on this scale
        sgs  <- .cov2Cor(sg)
        mus  <- x%*%alpha
        muw[,corColC] <- mus[,corColC]
        
        tmp <- .conditionalMVNRcpp(w, mus, sgs, cdex = cdex, S)
        mus    <- tmp$mu
        sgs    <- tmp$vr
        muw[,cdex] <- mus
        
        if(S1 == 1){
          w[,ddex] <- matrix(rnorm(n,mus,sqrt(sgs[1])))
        } else {
          w[,ddex] <- .rMVN(n,mus,sgs)
        }
      } 
      yg <- w
      if(length(ccols) > 0){
        mmm <- yg[,ccols]
        mmm[mmm < 0] <- 0
        yg[,ccols]   <- mmm
      } 

    for(k in allTypes){    # predicting from w (not from yg)
      
      wk  <- which(typeCols == k)
      nk  <- length(wk)
      wo  <- which(wk %in% notOther)
      wu  <- which(typeCols[notOther] == k)
      wp  <- w[, wk, drop=F]
      
      groups <- NULL
      
      if( typeFull[wk[1]] == 'countComp' ){
        
        groups <- CCgroups[wk]
        nkk    <- max(groups)
        
        for(j in 1:nkk){
          
          wjk <- which(typeCols[wk] == k & CCgroups[wk] == j)
     #     wss <- which(wjk %in% notOther)
          wno <- which(wk %in% notOther)
          woo <- which(wk %in% other)
          www <- w[,wk]
          www[www < 0] <- 0
          
          www <- .gjamCompW2Y(www,notOther=wno)$ww
          
          if(COND){
            www <- sweep(www,1,CCsums[[j]],'*')
          } else {
            www <- sweep(www,1,ysum,'*')
          }
          yg[,wk] <- www
        }
        
      } else {
        
        if(typeFull[wk[1]] == 'fracComp') groups <- FCgroups[wk]
        
        tmp <- .gjamWLoopTypes(wo, type = typeFull[wk[1]], yy = yg[,wk,drop=F], 
                               wp, yp = yg[,wk,drop=F], cutg, 
                               censor, censorCA, censorDA, effMat[,wk,drop=F], groups, 
                               k, typeCols, notOther, wk = wk )
        yg[,wk] <- tmp[[2]] #[,wk]
        yg[,wk] <- .censorValues(censor,yg,yg)[,wk]
      }
    }
  }
    
    if(length(ccols) > 0){
      mmm <- muw[,ccols]
      mmm[mmm < 0] <- 0
      muw[,ccols] <- mmm
    }
    
  #  if(!is.null(effMat))yg <- yg*effMat
    
    yg[,condCols] <- ydataCond
    
    prPresent[yg > 0] <- prPresent[yg > 0] + 1
    
    ig <- ig + 1
    setTxtProgressBar(pbar,ig)
    
    ypred  <- ypred + yg
    ypred2 <- ypred2 + yg^2
    wcred  <- wcred + muw
    wcred2 <- wcred2 + muw^2
    
    
    ess[notOther,notOther]  <- .cov2Cor( t(agg)%*%covE%*%agg ) 
    emat[notOther,notOther] <- emat[notOther,notOther] + ess[notOther,notOther]
    
    if(FULL){
      ygibbs[kg,] <- as.vector(yg)
      wgibbs[kg,] <- as.vector(muw)
    }
    
    kg <- kg + 1
  }
  
  prPresent <- prPresent/nsim
  
  
  ematrix  <- emat/nsim
  
  
  yMu  <- ypred/nsim
  res  <- ypred2/(nsim - 1) - yMu^2
  res[res < tiny] <- tiny
  yPe <- sqrt(res) 
  
  wMu  <- wcred/nsim
  res  <- wcred2/(nsim - 1) - wMu^2
  res[res < tiny] <- tiny
  wSe <- sqrt(res)
  
  colnames(yMu) <- colnames(yPe) <- colnames(wMu) <- 
    colnames(wSe) <- ynames
  
  sdList <- list( yMu = yMu, yPe = yPe, wMu = wMu, wSe = wSe )
  
  piList <- NULL
  if(FULL){
    wLo <- matrix( apply(wgibbs,2,quantile,.05), n, S )
    wHi <- matrix( apply(wgibbs,2,quantile,.95), n, S )
    yLo <- matrix( apply(ygibbs,2,quantile,.05), n, S )
    yHi <- matrix( apply(ygibbs,2,quantile,.95), n, S )
    
    colnames(wLo) <- colnames(wHi) <- colnames(yLo) <- 
      colnames(yHi) <-  ynames
    piList <- list( wLo = wLo, wHi = wHi, yLo = yLo, yHi = yHi )
  }
    
  if(PLOT){
    
    oma <- c(0,0,0,0)
    mar <- c(4,4,2,1)
    tcl <- -0.5
    mgp <- c(3,1,0)
    
    par(oma = oma, mar = mar, tcl = tcl, mgp = mgp, bty='n')
    
    wy <- which(colnames(y) %in% y2plot & c(1:S) %in% notOther)
    t2plot <- typeNames[wy]
    allTypes <- unique(t2plot)
    
    mfrow <- .getPlotLayout(length(allTypes) + 1)
    par(mfrow=mfrow, bty='n', mar=c(1,2,3,1) )
    
    k   <- 0
    add <- F
    
    for(j in 1:length(allTypes)){
      
      wk <- which(typeNames == allTypes[j] & c(1:S) %in% notOther)
      ws <- colnames(y)[wk]
      wm <- which(colnames(yMu) %in% colnames(y)[wk])
      wk <- match(colnames(yMu)[wm],colnames(y))
      
      y1 <- y[,wk]
      if(min(y1) == max(y1))next
      y2 <- yMu[,wm]
      
      tmp <- .gjamPlotPars(type=allTypes[j],y1,y2)
      y1 <- tmp$y1; yp <- tmp$yp; nbin <- tmp$nbin; nPerBin <- tmp$nPerBin
      vlines <- tmp$vlines; xlimit <- tmp$xlimit; ylimit <- tmp$ylimit
      breaks <- tmp$breaks; wide <- tmp$wide; LOG <- tmp$LOG; POINTS <- F
      MEDIAN <- tmp$MEDIAN
      
      log <- ''
      if(LOG)log <- 'xy'
      
      if(LOG){
        wn <- which(y1 > 0 & yp > 0)
        y1 <- y1[wn]
        yp <- yp[wn]
      }
      
      tmp <- .bins4data(y1,nPerBin=nPerBin,breaks=breaks,LOG=LOG)
      breaks <- tmp$breaks
      bins   <- tmp$bins
      nbin   <- tmp$nbin
      
      if( !allTypes[j] %in% c('PA','CAT') ){
        ncc   <- max( c(100,max(y1)/20) )
        xy <- .gjamBaselineHist(y1,bins=bins,nclass=ncc)
        xy[2,] <- ylimit[1] + .8*xy[2,]*diff(ylimit)/max(xy[2,])
        plot(xy[1,],xy[2,],col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
             xlab='Observed',ylab='Predicted')
        polygon(xy[1,],xy[2,],border='tan',col='wheat')
        
      } else {
        y11 <- mean(y1)
        y00 <- 1 - y11
        x11 <- c(-.07,-.07,.07,.07,.93,.93,1.07,1.07,-.07)
        y11 <- c(0,y00,y00,0,0,y11,y11,0,0)
        plot(x11,y11,col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
             xlab='Observed',ylab='Predicted')
        polygon(x11,y11,border='tan',col='wheat')
      }
      abline(0,1,lty=2,lwd=3,col='brown')
      abline(h = mean(y1),lty=2,lwd=3,col='tan')
      
      add <- T
      
      tmp <- .plotObsPred(y1,yp,xlabel='Observed',ylabel='Predicted',nbin=nbin, 
                          nPerBin=nPerBin, xlimit=xlimit,ylimit=ylimit,
                          breaks=breaks, wide=wide, LOG=LOG, fill='lightblue', 
                          box.col='darkblue',POINTS=F, MEDIAN=MEDIAN, add=add)
      if(length(vlines) > 0)abline(v=vlines,lty=2)
      
      tt <- allTypes[j]
      if(length(ws) == 1)tt <- paste(ws,tt,sep='-')
      
      lab <- paste(letters[j],') ',tt, sep='')
      .plotLabel( lab,above=T )
    }
    
    yp  <- colMeans(yMu)
    wy  <- match(colnames(yMu),colnames(y))
    
    .plotObsPred( colMeans(y[,wy]),yp, xlabel='Observed',
                  xlimit=NULL, ylimit=NULL,
                  breaks=breaks, wide=wide, LOG=LOG, fill='lightblue', 
                  box.col='darkblue', POINTS=T, ptcol='darkblue')
    abline(0, 1,lty=2,lwd=3,col='brown')
    abline(h = mean(y1),lty=2,lwd=3,col='tan')
    .plotLabel( paste(letters[j+1],') By Species',sep=''),above=T )
    
  }
    
  list( x = x, sdList = sdList, piList = piList, prPresent = prPresent,
        ematrix = ematrix)
}

   
.gjamUpdateBetaPrior <- function(WIX,IXX,sg,alpha,loBeta,hiBeta,...){
  
  S <- ncol(WIX)
  Q <- nrow(WIX)
  smat <- kronecker(sg,IXX)
  
  tmp <- .tnormMVNmatrixRcpp(avec = matrix(alpha,1), muvec = matrix(WIX,1), 
                             smat = smat, 
                             lo = matrix(loBeta,1), hi = matrix(hiBeta,1))
  tmp <- matrix(tmp,nrow(alpha),ncol(alpha))
  
  tmp[!is.finite(tmp)] <- alpha[!is.finite(tmp)]
  list(bg = tmp)
}
.gjamUpdateTheta <- function(w,tg,cutLo,cutHi,ordCols,holdoutN,
                             holdoutIndex,minOrd,maxOrd){
  
  word <- w[,ordCols,drop=F]
  ncut <- ncol(tg)
  nc   <- ncut - 1
  n    <- nrow(w)
  nk   <- length(ordCols)
  
  c1 <- cutLo[,1]
  c2 <- cutLo[,2]
  c3 <- cutHi[,1]
  c4 <- cutHi[,2]
  
  if(holdoutN > 0){
    word <- word[-holdoutIndex,]
    ss   <- seq(0,(nk-1)*n,by=n)
    wh <- as.vector( outer(holdoutIndex,ss,'+') )
    c1 <- c1[-wh]
    c2 <- c2[-wh]
    c3 <- c3[-wh]
    c4 <- c4[-wh]
  }
  
  cmin <- .byRcpp(as.vector(word),c1,c2,fun='min')
  cmax <- .byRcpp(as.vector(word),c1,c2,fun='max')
  
  cmin[!is.finite(cmin[,1]),1] <- -10
  cmin[,2] <- 0
  cmax[,1] <- 0
  cmax[cmax == -Inf] <- Inf
  
  tmp <- .interpRows(cmax,startIndex=minOrd+1,endIndex=maxOrd-1,
             INCREASING=T,minVal=0,maxVal=Inf,
             defaultValue=NULL,tinySlope=.001)
  
  cmax[!is.finite(cmax)] <- tmp[!is.finite(cmax)]
  
  ww <- which(!is.finite(cmin) & is.finite(cmax),arr.ind=T)
  if(length(ww) > 0){
    w0 <- ww
    w0[,2] <- w0[,2] - 1
    cmin[ww] <- runif(nrow(ww),cmax[w0],cmax[ww])
  }
  
 # ww <- which(is.finite(cmax))
  
  clo <- cmax[drop=F,,-nc]
  chi <- cmin[drop=F,,-1]
  clo[,1] <- -1
  
  ww <- which(is.finite(clo))
  cl <- clo[ww]
  ch <- chi[ww]
  wc <- which(cl > ch,arr.ind=T)
#  if(length(wc) > 0)stop(print(wc))
  cl[cl > ch] <- ch[cl > ch]
  
  chi[ww] <- .tnorm(length(ww),cl,ch,cl,3)
  chi[,1] <- 0
  cmax <- cbind(-Inf,chi,Inf)

  
  cmax[,ncut] <- Inf
  if( ncol(cmax) > max(maxOrd) )cmax[ cbind(1:nk,maxOrd+1) ] <- Inf
  
  wmin <- which(minOrd > 1)
  if(length(wmin) > 0){
    for(j in wmin)cmax[j,2:c(minOrd[j]+1)] <- 0:(minOrd[j] - 1)
  }
  cmax
}

.censorValues <- function(censor,y,yp){
  
  mm  <- length(censor)
  if(mm == 0)return(yp)
  
  if(mm > 0){
    for(m in 1:mm){
      wc  <- censor[[m]]$columns
      nc  <- ncol( censor[[m]]$partition )
      ym  <- yp[,wc,drop=F]
      cp  <- censor[[m]]$partition
      for(k in 1:nc){
  #      kcens <- which(y[,wc] == cp[1,k])
        wlo <- which( ym >  cp[2,k] & ym <  cp[3,k])
        ym[wlo] <- cp[1,k]
      }
      yp[,wc] <- ym
    }
  }
  yp
}


.gjamWLoopTypes <- function(wo, type, yy, wp, yp, cutg, censor, 
                            censorCA, censorDA,
                            effMat, groups, k, typeCols, notOther, wk = wo ){
  
  if( type == 'continuous')return( list(yy,yp) )   # w = y
  
  nk  <- ncol(wp)
  wkk <- c(1:nk)
  n  <- nrow(wp)
  
  if( type == 'ordinal' ){
    for(s in 1:nk)yp[,s] <- findInterval(yp[,s],cutg[s,]) - 1
    return( list(wp,yp) )
  }
  
  if( type == 'presenceAbsence' ){
    yp[yp > 0]  <- 1 
    yp[yp <= 0] <- 0
    return( list(wp,yp) )
  }
  
  if( type == 'contAbun' ){
    yp[yp < 0]    <- 0
    return( list(wp,yp) )
  }
  
  if( type == 'discAbun' ){
    yp[yp < 0] <- 0
    if(length(censorDA) > 0) wp[-censorDA] <- yy[-censorDA]
  #  if(!is.null(effort)){
  #    if(!is.matrix(effort$values)){
        yp <- yp*effMat
  #    } else {
  #      we  <- which(effort$columns %in% wk)
  #      wf  <- match(effort$columns[we],wk)
  #      yp[,wf] <- yp[,wf]*effort$values[,we]
  #    }
  #  }
    return( list(wp,yp) )
  }

  if( type == 'categorical' ){ ## only prediction
    
    ntt <- max( groups )
    
    for(i in 1:ntt){  
      
      if(ntt == 1){
        wki <- wkk
      } else {
        wki <- which( groups == i )
      }
      nki  <- length(wki)
      wko  <- wki
      wmax <- apply( yp[,wko],1, which.max) 
      
      yp[,wki] <- 0
      yp[,wki][ cbind(1:n,wmax) ] <- 1
    }
    return( list(wp,yp) )
  }
  
  if( type == 'countComp' ){  ##  w and y 
    
    ntt <- max( groups )
    ww  <- wp                 
    ww[ww < 0] <- 0     
    yp[yp < 0] <- 0
    
    for(i in 1:ntt){  ## normalize w and y 
      
      if(ntt == 1){
        wki <- wkk
      } else {
        wki <- which( groups == i )
      }
      
      io <- which(wki %in% wo)
      wc <- .gjamCompW2Y(ww[,wki], notOther=io)$ww
      wp[,wki][wp[,wki] > 0] <- wc[wp[,wki] > 0]
      
      yp[,wki] <- .gjamCompW2Y(yp[,wki],notOther=io)$ww
      ysum     <- rowSums(yy[,wki])
      yp[,wki] <- round( sweep(yp[,wki],1,ysum,'*'), 0) 
    }
    return( list(wp,yp) )
  }
  
  ## fracComp: w and y 
  
  ntt <- max( groups )
  
  wy     <- which(yy > 0)  
  wp[wy] <- yy[wy]     
  yp[yp < 0] <- 0
  
  for(i in 1:ntt){  ## normalize w and y 
    
    if(ntt == 1){
      wki <- wkk
    } else {
      wki <- which(groups == i)
    }
    
    io <- which(wki %in% wo)
    yp[,wki] <- .gjamCompW2Y(yp[,wki],notOther=io)$ww
  }
  return( list(wp,yp) )
}
  


.gjamWcatLoop <- function(y, ws, mus, sgs, notOther, plo, phi, groups, 
                          REDUCT = F){
  
  # if REDUCT, sgs is length-S sigvec
  # if !REDUCT, sgs is css[notOther,notOther]
  
  ntt <- max( groups )
  n <- nrow(y)
  
  for(i in 1:ntt){  
    
    wki <- which(groups == i)
    nki <- length(wki)
    wko <- wki[wki %in% notOther]
    
    w0 <- apply( ws[,wko]*(1 - y[,wko]),1, max) # max(w, 0) for y = 0
    w1 <- apply( ws[,wko]*y[,wko],1, max)       # w for y = 1
    w0[w0 < 0] <- 0                             # when y is reference
    
    si <- sample(wko)
    
    for(s in si){
      
      y1         <- which(y[,s] == 1)
      plo[-y1,s] <- -500
      phi[y1,s]  <- 500
      plo[y1,s]  <- w0[y1]
      phi[-y1,s] <- w1[-y1]
      
      if(REDUCT){
        
        ws[,s] <- .tnorm(n,plo[,s],phi[,s],mus[,s],sqrt(sgs[s]))
        
      } else {
        sm  <- which(notOther == s)
        tmp <- .conditionalMVNRcpp(ws[,notOther], mus[,notOther], 
                                   sgs, sm)
        mue <- tmp$mu
        vr  <- max(tmp$vr,1e-8)
        ws[,s] <- .tnorm(n,plo[,s],phi[,s],mue,sqrt(vr))
      }
      
      w1[y1]  <- ws[y1,s]    #new w for y = 1
      w0[-y1] <- apply( ws[-y1,wki]*(1 - y[-y1,wki]),1, max)
      
  #    print( round( cbind(w0,w1,ws[,wko])[1:10,], 3) )
    }
  }
  list(w = ws, plo = plo, phi = phi)
}

.gjamWcatLoop2 <- function(y, ws, mus, sgs, notOther, plo, phi, groups, 
                           REDUCT = F){
  
  # if REDUCT, sgs is length-S sigvec
  # if !REDUCT, sgs is css[notOther,notOther]
  
  ntt <- max( groups )
  n <- nrow(y)
  
  for(i in 1:ntt){  
    
    wki <- which(groups == i)
    nki <- length(wki)
    wko <- wki[wki %in% notOther]
    
    # should eliminate after first iteration, use plo, phi:
  #  w0 <- apply( ws[,wko]*(1 - y[,wko]),1, max ) # max(w, 0) for y = 0
    w1 <- apply( ws[,wko]*y[,wko],1, max)        # w for y = 1
 #   w0[w0 < 0] <- 0                              # when y is reference
    
    so <- match(wko,notOther)
    
    for(s in wko){   
      
      y1 <- which(y[,s] == 1)
      
      #   if(length(y1) == 0)next
      sm  <- which(notOther == s)   #index in sgs = sg[notOther,notOther]
      sn  <- so[so != sm]           #index in sgs for so
      qs  <- wko[wko != s]          
  #    plo[y1,s]  <- w0[y1]
      
      if(REDUCT){
        ws[y1,s] <- .tnorm(length(y1),plo[y1,s],phi[y1,s],
                           mus[y1,s],sqrt(sgs[s]))
      } else {
        tmp <- .conditionalMVNRcpp(ws[y1,notOther], mus[y1,notOther], sgs, sm)
        mue <- tmp$mu
        vr  <- max(tmp$vr,1e-8)
        ws[y1,s] <- .tnorm(length(y1),plo[y1,s],phi[y1,s],mue,sqrt(vr))
      }
      
      w1[y1] <- ws[y1,s]        # w for y = 1
      phi[y1,wki] <- w1[y1]
      phi[y1,s]   <- 500
      
      if(REDUCT){     # the zeros
        tmp <- .tnorm(length(y1)*length(qs),plo[y1,qs],phi[y1,qs],
                      mus[y1,qs],sqrt(sgs[s]))
      } else {
        tmp <- .tnormMVNmatrixRcpp(ws[y1,notOther],mus[y1,notOther],
                                   smat=sgs, plo[y1,notOther], 
                                   hi=phi[y1,notOther],
                                   whichSample=so)[,sn,drop=F]
      }
      ws[y1,qs] <- tmp
      
      ###########
      if(length(sn) > 0)tmp <- apply( tmp, 1, max ) #########
      tmp[tmp < 0] <- 0
      plo[y1,s]  <- tmp
    }
    ##############
    s <- wki[!wki %in% wko]   #  y = 1 is ref class
    y1 <- which(y[,s] == 1)
    tmp <- .tnormMVNmatrixRcpp(ws[y1,notOther],mus[y1,notOther],
                               smat=sgs, plo[y1,notOther], 
                               hi=phi[y1,notOther],
                               whichSample=so)
    ws[y1,wko] <- tmp[,so]
    #############
  }
  list(w = ws, plo = plo, phi = phi)
}



.gjamWLoop <- function(ws,mus,sgs,wkk,lo,hi,sampW, indexW,
                       byCol=T, byRow=F){
  
  n <- nrow(lo)
  tiny <- .00001
  
  if(byCol){
    
    iss <- wkk[wkk %in% indexW]
    
    for(s in iss){
      
      rs <- which(sampW[,s] > 0)
      ls <- lo[drop=F,rs,s]
      hs <- hi[drop=F,rs,s]
      
      tmp <- .conditionalMVNRcpp(ws[drop=F,rs,],mus[drop=F,rs,],sgs,s)
      mu  <- tmp$mu
      vr  <- max(tmp$vr,tiny)
      tmp <- .tnorm(length(rs),ls,hs,mu,sqrt(vr))
      
      wl  <- which(tmp == ls)
      if(length(wl) > 0) tmp[wl] <- ls[wl] + tiny*(ls[wl])
      
      wl  <- which(tmp == hs)
      if(length(wl) > 0) tmp[wl] <- hs[wl] - (1 - tiny)*hs[wl]
      
      ws[rs,s] <- tmp
    }
    return(ws)
  }
  
  for(i in indexW){
    
    rs  <- which(sampW[i,] > 0)
    rs  <- rs[rs %in% wkk]
    ws[i,rs] <- .tnormMVNmatrixRcpp(ws[drop=F,i,], mus[drop=F,i,],
                               smat=sgs, lo[drop=F,i,], hi[drop=F,i,],
                               whichSample=rs)[,rs]
  }
  ws
}


.setContrasts <- function(xx){
  
  # contrasts where each level is compared to the reference level
  # data must have an attribute for 'reference' class assigned as, e.g.,
  # attr(xdata$soil,'reference') <- 'reference'
  # where xx is xdata$soil and 'reference' is the name that appears in xx
  
  levs  <- attr(xx,'levels') 
  nl    <- length(levs)
  ml    <- nl - 1
  ref <- levs[1]
  
 # if(is.null(ref))return( list(levs = as.character(levs), cont = contrasts(xx)) )
  
  intType <- attr(xx,'intType')
  
  if(is.null(intType))intType <- 'ref'
    
  wr  <- which(levs == ref)
  
  cj <- matrix(-1/nl,ml,ml)
  diag(cj) <- ml/nl
  rownames(cj) <- levs[-wr]
  colnames(cj) <- levs[-wr]
  
  rj <- rep(-1/nl,ml)
  cj <- rbind(rj,cj)
  rownames(cj)[1] <- ref
  
  levs <- as.character(levs)
  
  cj <- cj[drop=F,levs,]
  if(intType == 'ref'){
    cj[cj > 0] <- 1
    cj[cj < 0] <- 0
  }
  list(levs = levs, cont = cj)
}



.gjamXY <- function(formula, xdata, y, typeNames, notStandard, 
                    checkX=T, xscale = NULL){
  
  n        <- nrow(xdata)
  S        <- ncol(y)
  snames   <- colnames(y)
  facNames <- character(0)
  
  original <- colnames(xdata)
  
  # colnames(xdata) <- .replaceString(colnames(xdata),'-','')
  # colnames(xdata) <- .replaceString(colnames(xdata),'_','')
  # colnames(xdata) <- .replaceString(colnames(xdata),' ','')
  
  # new <- colnames(xdata)
  # xdataNames <- rbind(original, new)
  xdataNames <- original
  
  if(!is.null(notStandard)){
    notStandard <- .replaceString(notStandard,'-','')
    notStandard <- .replaceString(notStandard,'_','')
    notStandard <- .replaceString(notStandard,' ','')
  }
  
  form <- attr( terms(formula), 'term.labels' )
  form <- .replaceString(form,'-','')
  form <- .replaceString(form,'_','')
  form <- paste0(form,collapse=' + ')
  
  formula <- as.formula( paste('~',form) )
  
  t1 <- attr( terms(formula), 'term.labels' )
  t2 <- sapply(xdata,is.factor)
  t2 <- t2[names(t2) %in% t1]
  
  facNames <- names(t2)[t2]   #facNames
  xxn      <- names(t2)[!t2]
  
  xmean <- xsd <- NULL
  
  xnames <- colnames(xdata)
  standX <- xxn  # 
  
  if(is.null(xscale) & length(standX) > 0){
    sc    <- which( xnames %in% standX)
    sc    <- match(xnames[sc],colnames(xdata))
    xmean <- colMeans(xdata[,sc,drop=F],na.rm=T)
    xsd   <- apply(xdata[,sc,drop=F],2,sd,na.rm=T)
    xscale <- rbind(xmean,xsd)
    xscale['xmean',abs(xscale['xmean',]) < 1e-15] <- 0
  }
  
  
  
  standX <- standX[!xxn %in% notStandard]
  
  
  if(length(standX) > 0){
    sc    <- which( xnames %in% standX )
    ss <- xnames[sc]
    xdata[,ss] <- (xdata[,ss] - matrix(xscale['xmean',ss],n,length(ss),byrow=T))/
      matrix(xscale['xsd',ss],n,length(ss),byrow=T)
    xmean <- xscale['xmean',]
    xsd   <- xscale['xsd',]
  }
  
  factorList <- contrast <- vector('list',length = length(facNames))
  names(factorList) <- facNames
  
  if(length(facNames) > 0){
    
    rr <- 1
    
    for(j in 1:length(facNames)){
      
      wj <- which(names(xdata) == facNames[j])
      xf <- as.character(xdata[[wj]])
      
      wk <- which(xf == facNames[j])         # levels same name as factor
      if(length(wk) > 0){
        xf[wk] <- paste('ref',rr,sep='')
        xf     <- as.factor(xf)
        xdata[[wj]] <- xf
        rr <- rr + 1
      }
      
      cj <- attr(xdata[[wj]],'contrasts')
      contrast[[j]] <- cj
      tt <- .setContrasts(xdata[[wj]])$cont
      factorList[[j]] <- paste(facNames[j],colnames(tt),sep='')
      
      if(!is.null(cj))next                       # contrasts previously set

      contrast[[j]] <- tt
      attr(xdata[[wj]],'contrasts') <- tt
    }
    names(contrast) <- facNames
  }  
  
  tmp <- model.frame(formula,data=xdata,na.action=NULL)
  x   <- model.matrix(formula, data=tmp)
  
  colnames(x)[1] <- 'intercept'
  
  xnames <- colnames(x)
  snames <- colnames(y)
  Q      <- ncol(x)
  predXcols <- 2:Q
  
  isFactor <- character(0)
  if(length(facNames) > 0){
    
    for(j in 1:length(facNames)){
      ij <- grep(facNames[j],colnames(x))
      ij <- xnames[ij]
      ix <- grep(':',ij)
      if(length(ix) > 0)ij <- ij[-ix]
      isFactor <- c(isFactor,ij)
      print(paste('observations in factor',facNames[j]))
      print(colSums(x, na.rm=T)[ij])
    }
  }
  
  VIF <- isNonLinX <- designTable <- NULL
  isInt <- intMat <- numeric(0)
  
  # check design
  
  if(checkX & length(xxn) > 0){
    checkInt <- range(x[,1])
    if(checkInt[1] != 1 | checkInt[2] != 1)
      stop( paste('x[,1] must be intercept (ones)') )
    
    tmp <- .checkDesign(x[,c('intercept',xxn)])
    if(tmp$rank < tmp$p)stop( 'x not full rank' )
    VIF         <- tmp$VIF
    designTable <- tmp$designTable$table
  }
  
  if(Q > 2 & length(xxn) > 0){
    
    wx <- grep('^2',colnames(x),fixed=T)
    if(length(wx) > 0){
      mm <- unique(unlist(strsplit(colnames(x)[wx],'^2)',fixed=T)))
      mm <- .replaceString(mm,'I(','')
      mm <- match(mm,colnames(x))
      mat <- cbind(wx,mm,mm)
      colnames(mat) <- c('int','main1','main2')
      intMat <- mat
      isInt <- wx
      isNonLinX <- sort(unique( c(isNonLinX,mm,isInt)))
    }
    
    wx <- grep(':',colnames(x))
    if(length(wx) > 0){
      mm  <- matrix(unlist(strsplit(colnames(x)[wx],':')),ncol=2,byrow=T)
      mat <- matrix( match(mm,colnames(x)), ncol=2)
      mat <- cbind(wx,mat)
      colnames(mat) <- c('int','main1','main2')
      wx <- c( which(colnames(x) %in% mm),wx )
      isInt <- sort(c(isInt,wx))
      intMat <- rbind(intMat,mat)
    }
    if(!is.null(isInt))isNonLinX <- sort(unique( c(isNonLinX,isInt)))
  }
  
  standMat <- matrix(1,Q,S)
  colnames(standMat) <- snames
  rownames(standMat) <- xnames
  standMu <- standMat - 1
  
  xss <- colnames(xscale)
  
  if(length(xss) > 0){
    standMu[xss,]  <-  xscale['xmean',xss]
    standMat[xss,] <-  xscale['xsd',xss]
  }
  
  # standardize in interactions
  
  if(length(intMat) > 0){
    
    for(j in 1:nrow(intMat)){
      im <- intMat[j,]
      s1 <- s2 <- 1
      if( xnames[im[2]] %in% colnames(xscale) )s1 <- xscale['xsd',xnames[im[2]]]
      if( xnames[im[3]] %in% colnames(xscale) )s2 <- xscale['xsd',xnames[im[3]]]
      
      standMat[im[1],] <- s1*s2
    }
  }
  
  standRows <- which(standMat[,1] != 1 | standMu[,1] != 0)
  standRows <- standRows[!names(standRows) %in% notStandard]
  
  colnames(y) <- .replaceString(colnames(y),':','x')  
  colnames(y) <- .replaceString(colnames(y),' ','')  
  
  # check composition
  
  tiny <- 1 + 1e-10
  
  if('FC' %in% typeNames){
    
    groups <- attr(typeNames,'FCgroups')
    
    if(is.null(groups)){
      groups <- rep(0,S)
      groups[typeNames == 'FC'] <- 1
      attr(typeNames,'FCgroups') <- groups
    }
    
    ngg    <- max(groups)
    for(kk in 1:ngg){
      wf <- which(groups == kk)
      if(length(wf) == 0)stop( 'FC data must have > 1 column' )
      ww <- which(y[,wf] < 0)
      if(length(ww) > 0)stop( 'FC values cannot be < 0' )
      wr <- rowSums(y[,wf],na.rm=T)
      vv <- unique(wr)
      ww <- which(vv != 0 & vv > 1.01)
      if(length(ww) > 0){
        wx <- which(wr %in% vv)
        ii <- paste0(wx, collapse=', ')
        stop( paste('FC data must sum to zero (all absent) or one, check obs:',ii))
      }
    }
  }
  
  if('CC' %in% typeNames){
    wf <- which(typeNames == 'CC')
    if(length(wf) < 2)stop( 'CC data must have > 1 column' )
  }
  
  if(is.null(snames))snames <- paste('S',1:S,sep='-')
  if(is.null(xnames))xnames <- paste('x',1:Q,sep='-')
  
  snames <- sub('_','-',snames)
  xnames <- sub('_','-',xnames)
  
  colnames(y) <- snames
  colnames(x) <- xnames
  
  if(length(isNonLinX) == 0)isNonLinX <- NULL
  
  list(x = x, y = y, snames = snames, xnames = xnames, predXcols = predXcols,
       isInt = isInt, intMat = intMat, facNames = facNames,
       xdata = xdata,
       factorList = factorList, isFactor = isFactor, isNonLinX = isNonLinX, 
       designTable = designTable, xmean = xmean, xscale = xscale,
       standMu = standMu, standMat = standMat, standRows = standRows,
       contrast = contrast, xdataNames = xdataNames, formula = formula)
}


.gjamCompW2Y <- function(ww,notOther=c(1:(ncol(ww)-1))){
  
  pg <- .995
 # pgPrior <- c(9000, 278)
  
  n  <- nrow(ww)
  W  <- rowSums(ww[,notOther])
  wh <- which(W > pg)
  other <- c(1:ncol(ww))[-notOther]
  
  if(length(wh) > 0){
    contract <- (1 - (1 - pg)^(W[wh]/pg))/W[wh]
    ww[wh,]  <- ww[wh,]*contract        
  }
  
  ww[,other] <- 1 - rowSums(ww[,notOther])
  
#  pg <- rbeta(1,pgPrior[1] + n - length(wh),pgPrior[2] + length(wh) )
  
  list(pg = pg, ww = ww )
}

.imputX_MVN <- function(xx,yy,beta,xmiss,sinv,xprior=0,xbound=NULL,priorWT=1){
  
  # priorWT is inverse of variance
  
  wcol <- unique(xmiss[,2])
  S    <- nrow(sinv)
  Q    <- nrow(beta)
  
  if(is.null(xbound))xbound <- apply(xx,2,range,na.rm=T)
  
  for(j in wcol){
    
    wj <- xmiss[drop=F,xmiss[,2] == j,]
 #   if(!is.matrix(wj))wj <- matrix(wj,1,2)
    wr <- wj[,1]
    xp <- xprior[xmiss[,2] == j]
    
    bj <- matrix(beta[j,],1)
    bn <- matrix(beta[-j,],Q - 1)
    
 #   xn <- matrix(xx[wr,-j],length(wr))
    xn <- xx[drop=F,wr,-j]
    z <- yy[drop=F,wr,] - xn%*%bn
    datwt <- bj%*%sinv%*%t(bj)
    V     <- 1/( datwt + priorWT*datwt )
    v     <- z %*%sinv%*%t(bj) + xp*priorWT
    
    xx[wj] <- .tnorm(length(wr),xbound[1,j],xbound[2,j],v%*%V,sqrt(V))
  }
  xx
}

.interp <- function(y,INCREASING=F,minVal=-Inf,maxVal=Inf,defaultValue=NULL,
                   tinySlope=NULL){  #interpolate vector x

  if(is.null(defaultValue))defaultValue <- NA

  tiny <- .00001
  if(!is.null(tinySlope))tiny <- tinySlope

  y[y < minVal] <- minVal
  y[y > maxVal] <- maxVal

  n  <- length(y)
  wi <- which(is.finite(y))

  if(length(wi) == 0)return(rep(defaultValue,n))
  if(length(wi) == 1)ss <- tiny

  xx  <- c(1:n)
  z  <- y

  if(wi[1] != 1) wi <- c(1,wi)
  if(max(wi) < n)wi <- c(wi,n)

  ss <- diff(z[wi])/diff(xx[wi])

  ss[is.na(ss)] <- 0

  if(length(ss) > 1){
    if(length(ss) > 2)ss[1] <- ss[2]
    ss[length(ss)] <- ss[length(ss)-1]
  }
  if(INCREASING)ss[ss < tiny] <- tiny

  if(is.na(y[1]))  z[1] <- z[wi[2]] - xx[wi[2]]*ss[1]
  if(z[1] < minVal)z[1] <- minVal
  if(z[1] > maxVal)z[1] <- maxVal

  for(k in 2:length(wi)){

     ki <- c(wi[k-1]:wi[k])
     yk <- z[wi[k-1]] + (xx[ki] - xx[wi[k-1]])*ss[k-1]
     yk[yk < minVal] <- minVal
     yk[yk > maxVal] <- maxVal
     z[ki] <- yk
  }
  z
}

.interpRows <- function(x,startIndex=rep(1,nrow(x)),endIndex=rep(ncol(x),nrow(x)),
                       INCREASING=F,minVal=-Inf,maxVal=Inf,
                       defaultValue=NULL,tinySlope=.001){  
  #interpolate rows of x subject to increasing

  nn  <- nrow(x)
  p  <- ncol(x)
  xx <- c(1:p)

  if(length(minVal) == 1)minVal <- rep(minVal,nn)
  if(length(maxVal) == 1)maxVal <- rep(maxVal,nn)

  ni   <- rep(NA,nn)
  flag <- numeric(0)

  z <- x

  for(i in 1:nn){
    if(startIndex[i] == endIndex[i]){
      z[i,-startIndex[i]] <- NA
      next
    }
    z[i,startIndex[i]:endIndex[i]] <- .interp(x[i,startIndex[i]:endIndex[i]],
                                             INCREASING,minVal[i],maxVal[i],
                                             defaultValue,tinySlope)
  }
  
  z
}

.invertSigma <- function(sigma,sigmaerror=NULL,otherpar=NULL, REDUCT){
  
  if(REDUCT){
    sinv <- .invWoodburryArma(sigmaerror, otherpar$Z[otherpar$K,])
  } else {
    sinv <- chol2inv(chol( sigma ))
  }
  sinv
}

.invMatZero <- function(sgibbs,nsim=2000,snames,knames,index=NULL,
                        COMPRESS = F, REDUCT=F,
                        sigErrGibbs = NULL, kgibbs = NULL,
                        otherpar = NULL, alpha = .95){   
  # return conditional independence
  # if COMPRESS, sgibbs is as.vector(lower.tri(Sigma,diag=T) )
  # alpha: prob that covariance/inverse is not zero 
  
  S <- length(snames)
  
  if(is.null(index))index <- c(1:nrow(sgibbs))
  simIndex   <- sample(index,nsim,replace=T)
  
  if(!REDUCT){
    
    if(COMPRESS){
      tmp <- .expandSigmaChains(snames, sgibbs, otherpar, 
                                simIndex, sigErrGibbs, kgibbs, 
                                REDUCT=REDUCT, CHAINSONLY=T)$chainList$schain
      sgibbs <- tmp
    }
    S1 <- sqrt(ncol(sgibbs))
  #  if(is.null(colnames(sgibbs))){
  #    cnames <- paste('S',c(1:S1),sep='')
  #    cnames <- outer( cnames,cnames,paste,sep='_')
  #    colnames(sgibbs) <- cnames
  #    knames <- cnames
  #  }
  } else {
    N  <- otherpar$N
    r  <- otherpar$r
    S1 <- S
    SS <- matrix(0,S1,S1)
  }
  
  SK     <- length(knames)
  sindex <- match(knames,snames)
  mm     <- matrix(0,SK,SK)
  rownames(mm) <- colnames(mm) <- knames
  hiSS <- loSS <- hiSI <- loSI <- mm
  
  
  for(j in simIndex){
    if(!REDUCT){
      ss <- matrix(sgibbs[j,],S1,S1) #[sindex,sindex]
      si <- chol2inv(chol( ss ) ) 
      
    } else {
      Z  <- matrix(sgibbs[j,],N,r)
      ss <- .expandSigma(sigErrGibbs[j], S1, Z = Z, kgibbs[j,], REDUCT = T)
      si <- .invWoodburryArma(sigErrGibbs[j], Z[kgibbs[j,],])
    }
    
    ss <- ss[sindex,sindex]
    si <- si[sindex,sindex]
    
    hiSS[ss > 0] <- hiSS[ss > 0] + 1/nsim
    loSS[ss < 0] <- loSS[ss < 0] + 1/nsim
    hiSI[si > 0] <- hiSI[si > 0] + 1/nsim
    loSI[si < 0] <- loSI[si < 0] + 1/nsim
  }
  
  loMar <- which(loSS > alpha)
  hiMar <- which(hiSS > alpha)
  inMar <- which(loSS < alpha & hiSS < alpha)   # not different from zero
  
  loCon <- which(loSI > alpha)
  hiCon <- which(hiSI > alpha)
  inCon <- which(loSI < alpha & hiSI < alpha)
  
  inMarMat <- which(loSS < alpha & hiSS < alpha,arr.ind=T)
  inConMat <- which(loSI < alpha & hiSI < alpha,arr.ind=T)
  
  list( inMarMat = inMarMat, inConMat = inConMat )
}

.joinCharVec <- function(charVec,sep=''){
  
  if(length(charVec) == 1) return(charVec)
  xx <- charVec[1]
  for(j in 2:length(charVec))xx <- paste(xx,charVec[j],sep=sep)
  xx
}

.mapSetup <- function(xlim,ylim,scale=NULL,widex=10.5,widey=6.5){  
  
  #scale is x per inch
  #new means not a new plot
  
  if(is.null(scale))scale <- 1

  px   <- diff(xlim)/scale
  py   <- diff(ylim)/scale
  
  if(px > widex){
    dx <- widex/px
    px <- widex
    py <- py*dx
  }
  if(py > widey){
    dx <- widey/py
    py <- widey
    px <- px*dx
  }
    
  par(pin=c(px,py))
  
  invisible( c(px,py) )

}

.sameByColumn <- function(mat,fraction=F){
  
  nc <- ncol(mat)
  
  sameMat <- matrix(0,nc,nc)
  
  for(j in 2:nc){
    for(k in 1:(j - 1)){
       wj <- which(mat[,j] == mat[,k])
       sameMat[j,k] <- length(wj)
    }
  }
  fraction <- sameMat/nrow(mat)
  fraction[upper.tri(fraction, diag=T)] <- NA
  fraction
}

.modalValuesInArray <- function(mat,idim = 1){
  
  # modal values for each row (idim = 1) or column (idim = 2)
  
  as.numeric( apply(mat,idim,
                    function(x) names(which.max(table(x)))) )
}

.multivarChainNames <- function(rowNames,colNames){
  as.vector( t(outer(colNames,rowNames,paste,sep='_')) )
}

.multivarEmat <- function(bchains,covx,snames,orderB, nsim=500, 
                          alpha = .95){
  
  # SB - columns in bchains, excludes other
  # alpha - prob that cor is not zero
  
  Q  <- ncol(covx)
  SM <- length(snames)
  SB <- length(orderB)
  
  nc <- ncol(bchains)
  
  onames <- snames[orderB]
  bn     <- matrix( unlist(strsplit(colnames(bchains),'_')),nc,ncol=2,byrow=T)
  
  bcols <- ccols <- numeric(0)
  ccols <- character(0)
  for(j in 1:SB){
    wj    <- which(onames[j] == bn[,1])
    if(length(wj) == 0)wj <- which(onames == bn[,2])
    if(length(wj) == 0)next
    bcols <- c(bcols,min( wj ) )
    ccols <- c(ccols,colnames(bchains)[wj])
  }
  
  cnames <- colnames(bchains[,bcols])
  cnames <- .replaceString(cnames,now='intercept_',new='')
  cnames <- .replaceString(cnames,now='_intercept',new='')
  
  ng    <- nrow(bchains)  
  index <- c(round(ng/10,0):ng)
  
  jj <- sample(index,nsim,replace=T)
 # rj <- matrix(NA,nsim,SB*SB)
  
  bm <- matrix(0,SB,SB)
  colnames(bm) <- rownames(bm) <- cnames
  lo <- hi <- lm <- hm <- bm
  
  
  pbar <- txtProgressBar(min=1,max=nsim,style=1)
  ig   <- 1
  
  for(j in jj){
    
    bb     <- matrix( bchains[j,ccols],Q,SB) #includes omitSpec, but not other
    ss     <- t(bb)%*%covx%*%bb
    rr     <- .cov2Cor(ss) 
 #   rj[j,] <- as.vector( rr )
    bm <- bm + rr
    
    ri <- ginv(rr)
    
    lo[ rr < 0 ] <- lo[ rr < 0 ] + 1/nsim
    hi[ rr > 0 ] <- hi[ rr > 0 ] + 1/nsim
    lm[ ri < 0 ] <- lm[ ri < 0 ] + 1/nsim  # neg values
    hm[ ri > 0 ] <- hm[ ri > 0 ] + 1/nsim  # pos values
    
    setTxtProgressBar(pbar,ig)
    ig <- ig + 1
    
  }
  
  bm <- bm/nsim
  
  whichZero <- which(lo < alpha & hi < alpha,arr.ind=T) #not different from zero
  whConZero <- which(lm < alpha & hm < alpha,arr.ind=T)
  
  list(bm = bm, whichZero = whichZero, whConZero = whConZero)
}

.rMVN <- function (nn, mu, sigma){
  
  # nn - no. samples from one mu vector or nrow(mu) for matrix
  
  if(!is.matrix(mu)) mu <- matrix(mu,1)
  if(length(mu) == 1)mu <- matrix(mu,1,nrow(sigma))
  if(ncol(mu) == 1)  mu <- t(mu)
  
  m <- ncol(sigma)
  
  if(ncol(mu) != m)stop('dimension mismatch mu, sigma')
  
  if(nn > 1 & nrow(mu) == 1)mu <- matrix(mu,nn,m,byrow=T)
  
  if(nn != nrow(mu))stop('sample size does not match mu')
  
  testv <- try(svd(sigma),T)
  
  if( inherits(testv,'try-error') ){
    ev <- eigen(sigma, symmetric = TRUE)
    testv <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  } else {
    testv <- t(testv$v %*% (t(testv$u) * sqrt(testv$d)))
  }
  
  p <- matrix(rnorm(nn * m), nn) %*% testv
  p + mu
}
.sigCols <- function(mat,alpha=.025){
  #columns in mat that are different from zero
  
  tmp <- apply(mat,2,quantile,c(alpha,1 - alpha))
  pos <- which(tmp[1,] > 0)
  neg <- which(tmp[2,] < 0)
  list(pos = pos, neg = neg)
}
.omitChainCol <-
function(cmat,omitCols){
  
  #omitCols - characterVector
  
  keep <- c(1:ncol(cmat))
  ocol <- numeric(0)
  for(j in 1:length(omitCols)){
    ocol <- c(ocol,grep(omitCols[j],colnames(cmat)))
  }
  if(length(ocol) > 0)keep <- keep[-ocol]
  list(keep = keep, omit = ocol)
}

.outFile <- function(outfolder=character(0),file){
  paste(outfolder,file,sep='/')
}

.plotChainDensity <- function(cMat,vnames=NULL,ncolPlot=2,xlim=NULL,
                              xlab=' ',ylab=' ', LOG=F,inColor=F){
  
  if(is.null(vnames)){
    mp <- 1
  } else {
    mp  <- length(vnames)
  }
  
  knames <- vnames
  if(is.null(vnames))knames <- 'all'
  
  for(k in 1:length(knames)){
    
    kk <- knames[k]
    if(is.null(vnames))kk <- NULL
    
    if(!is.null(kk)){
      gg <- grep(kk,colnames(cMat))
      if(length(gg) == 0)next
    }
    
    tmp <- .chains2density(cMat,varName=kk, cut=3)
    xt  <- tmp$x
    yt  <- tmp$y
    chainMat <- tmp$chainMat
    
    nn <- nrow(chainMat)
    
    xrange <- quantile(xt,c(.005,.995))
    if(is.null(xlim))xlim <- range(xt)
    if(!LOG) plot(10,10,xlim=xlim,ylim=c(0,1.8*max(yt)),xlab=xlab,
                  ylab=ylab,cex=.01)
    if(LOG)  plot(10,10,xlim=xlim,ylim=c(0,1.8*max(yt)),xlab=xlab,
                  ylab=ylab,cex=.01,log='x')
    
    .plotLabel(kk,above=T)
    
    j1 <- 1
    if(knames[1] == 'intercept')j1 <- 2
    
    for(j in j1:nrow(xt)){
      
      xj <- xt[j,]
      yj <- yt[j,]
      
      cj <- cumsum(yj)
      cj <- cj/max(cj)
      ci <- xj[ findInterval(c(.05,.95),cj) ]
      
      label <- rownames(xt)[j]
      
      wm <- which.max(yj)
      lines(xt[j,],yj,lwd=2,col='brown')
      lines(range(xj),c(0,0),col='brown',lwd=1)
      
      polygon( c(xj,rev(xj)), c(yj,yj*0), border=.getColor('brown',.1),
               col = .getColor('brown',.1) )
      
      if(ncol(chainMat) < 25)text(xt[j,wm],1.1*yt[j,wm],label,
                                  srt=55,pos=4,cex=1)
    }
  }
}
.plotLabel <- function(label,location='topleft',cex=1.3,font=1,
                       above=F,below=F,bg=NULL){
  
  if(above){
    adj <- 0
    if(location == 'topright')adj=1
    title(label,adj=adj, font.main = font, font.lab =font)
    return()
  }
  if(below){
    adj <- 0
    if(location == 'bottomright')adj=1
    mtext(label,side=1,adj=adj, outer=F,font.main = font, font.lab =font,cex=cex)
    return()
  }
    
  if(is.null(bg)){
    tmp <- legend(location,legend=' ',bty='n')
  } else {
    tmp <- legend(location,legend=label,bg=bg,border=bg,text.col=bg,bty='o')
  }
  
  xt <- tmp$rect$left # + tmp$rect$w
  yt <- tmp$text$y
  
  pos <- 4
  tmp <- grep('right',location)
  if(length(tmp) > 0)pos <- 2
  
  XX <- par()$xlog
  YY <- par()$ylog
  
  if(XX)xt <- 10^xt
  if(YY)yt <- 10^yt
  
  text(xt,yt,label,cex=cex,font=font,pos=pos)
  
}

.bins4data <- function(obs, nPerBin=NULL, breaks=NULL, nbin=NULL, LOG=F){
  
  if( is.null(breaks) ){
    
    if(is.null(nbin))nbin <- 20
    
    br   <- range(obs[is.finite(obs)],na.rm=T)
    bins <- seq(br[1],br[2],length=nbin)
    if(LOG){
      yy <- obs[obs > 0]
      oo <- min( yy,na.rm=T )
 
      nper <- length(yy)/nbin
      nbb <- nper/length(yy)
      nbb <- seq(0,1,length=100)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      bins <- 10^quantile(log10(yy),nbb,na.rm=T)
      bins <- sort(unique(bins))
      br[1] <- .5*oo
     
    }
    if( !is.null(nPerBin) ){
      nbb <- nPerBin/length(obs)
      nbb <- seq(0,1,by=nbb)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      bins <- quantile(obs,nbb,na.rm=T)
      bins <- sort(unique(bins))
      
      db <- diff(bins)
      qo <- quantile(obs,c(.1,.9),na.rm=T)
      wb <- which( db/diff(range(qo)) < .02)
      wb <- wb[wb != 1]
      if(length(wb) > 0)bins <- bins[-wb]
      
      nbin <- length(bins)
    }
  } else {
    bins <- breaks
    nbin <- length(bins)
  }
  
  list(breaks = breaks, bins = bins, nbin = nbin)
}

.plotObsPred <- function(obs,yMean,ySE=NULL,nbin=NULL,nPerBin=NULL,breaks=NULL,
                         LOG=F,xlimit=NULL,ylimit=NULL,xlabel='Observed',
                         ylabel='Predicted',
                         ptcol=NULL, boxPerc = .6826895, whiskerPerc = .95,
                         fill=NULL, add=F, box.col='black', wide=NULL, POINTS=T,
                         MEDIAN=T){
  
  aa <- (1 - boxPerc)/2
  boxQuant <- c(aa, 1 - aa )
  aa <- (1 - whiskerPerc)/2
  whiskerQuant <- c(aa, 1 - aa )
  
  if(is.null(ptcol)){
    ptcol <- 'black'
 #   if(!is.null(nbin))ptcol <- 'grey'
  }
  if(length(ptcol) == 1)ptcol <- rep(ptcol,length(obs))
  
  if(is.null(xlimit))xlimit <- quantile(obs[is.finite(obs)],c(.01,.99),na.rm=T)
  if(is.null(ylimit))ylimit <- range(yMean[is.finite(yMean)],na.rm=T)
  
  xxx <- obs
  yyy <- yMean
  
  if(LOG){
    if(is.null(xlimit))xlimit <- range( obs[obs > 0],na.rm=T )
    if(is.null(ylimit))ylimit <- range( yMean[yMean > 0],na.rm=T )
    if(xlimit[1] <= 0)xlimit[1] <- .001
  }
  
  if(!POINTS){
    xxx <- xlimit[1]
    yyy <- ylimit[1]
  }
  
  if(!add){
    if(is.null(ylimit)){
      if(!LOG)plot(xxx,yyy,col=ptcol,cex=.03,xlab=xlabel,ylab=ylabel)
      if(LOG) plot(xxx,yyy,col=ptcol,cex=.03,xlab=xlabel,ylab=ylabel,log='xy')
    }
    if(!is.null(ylimit)){
      if(!LOG)plot(xxx,yyy,col=ptcol,cex=.03,xlab=xlabel,ylab=ylabel,
                   xlim=xlimit,ylim=ylimit)
      if(LOG) plot(xxx,yyy,col=ptcol,cex=.03,xlab=xlabel,ylab=ylabel,
                   xlim=xlimit,log='xy',ylim=ylimit)
    }
  }
  
  if(POINTS)points(xxx,yyy,pch=16,col=.getColor(ptcol,.5), cex=.5)
  
  if(!is.null(ySE)){
    ylo <- yMean - 1.96*ySE
    yhi <- yMean + 1.96*ySE
    for(i in 1:length(obs))lines(c(obs[i],obs[i]),c(ylo[i],yhi[i]),
                                 col='grey',lwd=2)
  }
  
  tmp    <- .bins4data(obs,nPerBin=nPerBin,breaks=breaks,LOG=LOG)
  breaks <- tmp$breaks
  bins   <- tmp$bins
  nbin   <- tmp$nbin
  
  if(is.null(wide))wide <- diff(bins)/2.1
  if(length(wide) == 1)wide <- rep(wide,nbin)
  minmax <- par('usr')[1:2]
  dff    <- diff(minmax)
  if(!LOG)wide[wide > dff/5] <- dff/5
  
  maxx <- 0
  last <- F
  
  for(k in 1:(nbin-1)){
    
    mb <- bins[k+1]
    if(mb >= xlimit[2]){
      last <- T
      mb   <- Inf
    }
    ok <- which(obs >= bins[k] & obs < mb)
    if(length(ok) == 0)next
    qk <- which(is.finite(yMean) & obs >= bins[k] & obs <= mb)
    q  <- quantile(yMean[qk],c(.5,whiskerQuant[1],boxQuant[1],
                               boxQuant[2],whiskerQuant[2]),na.rm=T)
    if(LOG)q[q <= 0] <- ylimit[1]
    
    ym <- q[1]
    xx <- mean(bins[k:(k+1)])      # bounded by bins
    if(!LOG){
  #    if(!is.null(nPerBin) & MEDIAN)xx <- median(obs[ok],na.rm=T)
      if(MEDIAN)xx <- median(obs[ok],na.rm=T)
    } else {
      xx <-  sqrt( prod(bins[k:(k+1)]) )
    }
    points(xx,q[1],pch=3,col=box.col)
    yy    <- q[c(2,5)]
    yy[1] <- max(c(yy[1],ylimit[1]),na.rm=T) + .0000001
    yy[2] <- max(yy)
    
    yy1 <- q[3]
    yy1 <- max(yy1,ylimit[1],na.rm=T) + .00000001
    yy2 <- max(yy1,q[4])
    
    minx <- xx - .4*(xx - bins[k])
    maxx <- xx + .4*(mb - xx)
    
    figRange <- par('usr')
    totalwide <- (maxx - minx)/diff(figRange[1:2])
    
    if(is.null(nPerBin)){

      if(maxx >= xlimit[2])maxx <- xlimit[2]
      
      if(LOG){
        
        dx   <- log10(bins[k+1]) - log10(xx)
        maxx <- 10^(log10(xx) + .2*dx)
        if(k == 1){
          dx <- -log10(xlimit[1]) + log10(xx)
        } else {
          dx <- -log10(bins[k-1]) + log10(xx)
        }
        minx <- 10^(log10(xx) - .2*dx)
        if(minx < xlimit[1])minx <- xlimit[1]
        totalwide <- (log10(maxx) - log10(minx))/diff(figRange[1:2])
      }
    
      
      rect(minx,yy1,maxx,yy2,col=fill,border=box.col)
      lines(c(minx,maxx),c(ym,ym),lwd=2,col=box.col)
      xx <- mean(c(minx,maxx))
    }
    
    if(!is.null(nPerBin)){
      qo <- quantile(obs[ok],c(.3,.7,.25,.75),na.rm=T)
      if(qo[1] == qo[2] | !MEDIAN)qo <- c(xx-.2*wide[k],
                                          xx+.2*wide[k],xx-.3*wide[k],
                                          xx+.3*wide[k])
      rect(qo[1],yy1,qo[2],yy2,col=fill,border=box.col)
      lines(c(qo[3],qo[4]),c(ym,ym),lwd=2,col=box.col)
      lines(rep(mean(qo[1:2]),2),yy,lwd=2,col=box.col)
    } else {
      lines(c(xx,xx),yy,lwd=2,col=box.col)
    }
    if(last)break
  }
  
  invisible( bins )
}

.predictY2X_linear <- function(xpred,yy,bb,ss,sinv=NULL, 
                               priorIV = diag(1e-10,ncol(xpred)), 
                               priorX = matrix(0,ncol(xpred)), 
                               predCols = c(2:ncol(xpred)),REDUCT, lox, hix){
  
  #inverse prediction for multivariate linear in x
  
  prX <- priorX[predCols]
  if(!is.matrix(prX))prX <- matrix(prX)
  
  nn <- nrow(yy)
  notPred <- c(1:ncol(xpred))[-predCols]
  
  bp <- matrix(bb[drop=F,predCols,],length(predCols))
  
  if(length(notPred) > 0){
    bn <- matrix(bb[notPred,],length(notPred))
    yy <- yy - xpred[,notPred]%*%bn
  }
  pp <- length(predCols)
  
  if(is.null(sinv))sinv <- chol2inv(chol(ss))
  
  bs <- bp%*%sinv
  
  V <- chol2inv(chol( bs%*%t(bp) + priorIV[predCols,predCols] ) )
  v <- yy%*%t(bs) + matrix( priorIV[predCols,predCols] %*% prX,nn,pp,byrow=T)
  mu <- v%*%V
  
  qq <- ncol(mu)
  
  if(qq > 1){
    xpred[,predCols] <- .tnormMVNmatrixRcpp(avec=xpred[,predCols],muvec=mu,smat=V,
                               lo=matrix(lox[predCols],nn,qq,byrow=T),
                               hi=matrix(hix[predCols],nn,qq,byrow=T))
  } else {
    xpred[,predCols] <- .tnorm(nn,lox[predCols],hix[predCols], mu,sqrt(V))
  }
  xpred
}

.predictY2X_nonLinear <- function(xx,yy,bb,ss,priorIV = diag(1e-10,ncol(xx)), 
                                  priorX=matrix(0,ncol(xx)),predCols=c(2:ncol(xx)),
                                  isInt=NULL,intMat=NULL, isFactor=NULL,
                                  factorList=NULL, contrast = contrast, lox, hix){
  
  #inverse prediction for multivariate nonlinear in x and factors, metropolis
  
  iFcol  <- NULL
  priorX <- priorX[predCols]
  if(!is.matrix(priorX))priorX <- matrix(priorX)
  
  nn <- nrow(yy)
  intercept <- xx[,1]
  
  xnew <- xx
  
  xv <- as.vector(xx[,predCols])
  nv <- length(xv)
  lo <- rep(lox[predCols],each=nn)
  hi <- rep(hix[predCols],each=nn)
  xnew[,predCols] <- .tnorm(nv,lo,hi,xv,.01)
  
 # xnew[,predCols] <- .rMVN(nn,xx[,predCols],diag(.01,length(predCols)))
  
  if(length(isFactor) > 0){          # all factors, main effects
 #   xnew[,isFactor] <- 0
    np <- length(factorList)
    for(k in 1:np){
      nf <- length(factorList[[k]]) + 1
      tm <- contrast[[k]][sample(nf,nn,replace=T),]
      xnew[,factorList[[k]]] <- tm
    }
    iFcol <- match(isFactor,colnames(xx))
  }
  
  if(length(intMat) > 0){     # some of the nlin terms interactions?
  #  pindex <- unique(as.vector(intMat[,-1,drop=F]))
  #  if(length(iFcol) > 0)pindex <- pindex[!pindex %in% iFcol]  #  not factors
  #  if(length(pindex) > 1)xnew[,pindex] <- .rMVN(nn,xx[,pindex],diag(.01,length(pindex)))
  #  if(length(pindex) == 1)xnew[,pindex] <- rnorm(nn,xx[,pindex],sqrt(.01))
    xnew[,intMat[,1]] <- xnew[,intMat[,2]]*xnew[,intMat[,3]]
  }
  
  pnow <- .dMVN(yy,xx%*%bb,smat=ss,log=T)
  pnew <- .dMVN(yy,xnew%*%bb,smat=ss,log=T)
  
  a  <- exp(pnew - pnow)
 # wa <- which(a == 1)    # rounding error when pr approx 0
 # if(length(wa) > 0){
 #   a[a == 1 & pnew < pnow] <- 0
 # }
    
  z  <- runif(nn,0,1)
  wa <- which(z < a)
  xx[wa,] <- xnew[wa,]
  
  list(x = xx, accept = length(wa))
}
.predVsObs <- function(true,p,xlim=range(true),ylim=range(p,na.rm=T),xlab=' ',
                       ylab=' ', colors=rep(1,length(true)),lwd=2,add=F){ 
	
  #true  - length n vector of obs or true values
  #p - ng by n matrix of estimates
  
  if(!is.matrix(p))p <- matrix(p,ncol=1)
  
  nn <- length(true)
  y  <- apply(p,2,quantile,c(.5,.025,.975))
  
  if(!add)plot(true,y[1,],xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,col=colors,pch=3,lwd=lwd)
  points(true,y[1,],col=colors,pch=3,lwd=lwd)

  for(j in 1:nn)lines(c(true[j],true[j]),y[2:3,j],col=colors[j],lwd=lwd)
  abline(0,1,lty=2)
  
  invisible(y)
}
.processPars <- function(xgb,xtrue=numeric(0),CPLOT=F,DPLOT=F,
                         sigOnly = F,burnin=1,xlimits = NULL){  

  #xg      - matrix of gibbs chains
  #xtrue   - true values (simulated data)
  #CPLOT   - if T, plot chains
  #DPLOT   - if T, plot density
  #burnin  - analyze chains > burnin
  #xlimits - xlimits for plot
  #sigOnly - plot only parameters that 95% CI does not include 0
  
  if(!is.matrix(xgb))xgb <- matrix(xgb,ncol=1)
  if(is.null(colnames(xgb)))colnames(xgb) <- paste('V',c(1:ncol(xgb)),sep='-')
  
  NOPARS <- F
  
  if(sigOnly){
    wi   <- grep('intercept',colnames(xgb))      #extract covariates for plotting
    btmp <- xgb
    if(length(wi) > 0){
    	btmp <- xgb[,-wi]
      if(length(xtrue) > 0)xtrue <- xtrue[-wi]
    }

    wq   <- apply(btmp,2,quantile,c(.025,.975),na.rm=T)  #extract parameters != 0
    wq   <- which(wq[1,] < 0 & wq[2,] > 0)
    
    if(length(wq) == ncol(btmp))NOPARS <- T
    if(NOPARS) warning('no significant pars to plot')
    if(length(wq) > 0 & !NOPARS){
      xgb  <- btmp[,-wq]
      if(length(xtrue) > 0)xtrue <- xtrue[-wq]
    }
   }

  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  if(burnin > 1){
  	     if(burnin > (nrow(xgb) + 100))stop("burnin too large")
  	     xgb <- xgb[-c(1:burnin),]
  }
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  nc <- ncol(xgb)
  nf <- round(sqrt(nc),0)

  out <- t(rbind(apply(xgb,2,mean,na.rm=T),apply(xgb,2,sd,na.rm=T),
                 apply(xgb,2,quantile,c(.025,.975),na.rm=T)))
  if(!is.null(colnames(xgb)))rownames(out) <- colnames(xgb)
  colnames(out) <- c('estimate','se','0.025','0.975')
  if(length(xtrue) > 0){
    out <- cbind(out,xtrue)
    colnames(out) <- c('estimate','se','0.025','0.975','true value')
  }

  if(CPLOT | DPLOT)par(mfrow=c((nf+1),nf),mar=c(4,2,2,2))
  if(CPLOT & DPLOT)par(mfrow=c((nf+1),nc),mar=c(4,2,2,2))

  if(CPLOT & !NOPARS){
      for(j in 1:nc){
       plot(xgb[,j],type='l')
       abline(h=out[j,],lty=2)
       if(length(xtrue) > 0)abline(h=xtrue[j],col='red')
       abline(h = 0, col='grey',lwd=2)
       title(colnames(xgb)[j])
     }
  }
  xlims <- xlimits
  if(DPLOT & !NOPARS){
      for(j in 1:nc){
        xj <- density(xgb[,j])
        if(is.null(xlimits))xlims <- range(xj$x)
        plot(xj$x,xj$y,type='l',xlim=xlims)
        abline(v=out[j,],lty=2)
        if(length(xtrue) > 0)abline(v=xtrue[j],col='red')
        title(colnames(xgb)[j])
     }
  }
  list(summary = signif(out,4)
)

}
.replaceString <- function(xx,now='_',new=' '){  #replace now string in vector with new
  
  ww <- grep(now,xx,fixed=T)
  if(length(ww) == 0)return(xx)
  
  for(k in ww){
    s  <- unlist( strsplit(xx[k],now,fixed=T) )
    ss <- s[1]
    if(length(s) == 1)ss <- paste( ss,new,sep='')
    if(length(s) > 1)for(kk in 2:length(s)) ss <- paste( ss,s[kk],sep=new)
    xx[k] <- ss
  }
  xx
}

.joinCharVec <- function(charVec,sep=''){
  
  if(length(charVec) == 1) return(charVec)
  xx <- charVec[1]
  for(j in 2:length(charVec))xx <- paste(xx,charVec[j],sep=sep)
  xx
}

.buildYdata <- function(ydata, ytypes){
  
  # when y has factors, data.frame to matrix
  
  S <- ncol(ydata)
  
  wd <- which(duplicated(colnames(ydata)))
  if(length(wd) > 0){
    warning('duplicated colummn names in ydata')
    for(k in 1:length(wd)){
      dname <- colnames(ydata)[wd[k]]
      wk    <- which(colnames(ydata) == dname)
      colnames(ydata)[wk] <- paste(dname,1:length(wk),sep='')
    }
  }
  
  original <- colnames(ydata)
  
  colnames(ydata) <- .replaceString(colnames(ydata),'-','')
  colnames(ydata) <- .replaceString(colnames(ydata),'_','')
  colnames(ydata) <- .replaceString(colnames(ydata),' ','')
  
  new <- colnames(ydata)
  ydataNames <- rbind(original,new)
  
  CCgroups  <- attr(ytypes,'CCgroups')
  FCgroups  <- attr(ytypes,'FCgroups')
  CATgroups <- attr(ytypes,'CATgroups')
  
#  if( is.matrix(ydata) ){
#    if(is.null(CCgroups))CCgroups = rep(0,ncol(ydata))
#    if(is.null(FCgroups))FCgroups = rep(0,ncol(ydata))
    
#    if( 'FC' %in% ytypes & is.null(attr(ytypes,'FCgroups')) ){
#      groups <- rep(0,S)
#      groups[ytypes == 'FC']  <- 1
#      attr(ytypes,'FCgroups') <- groups
#    }
#    if( 'CC' %in% ytypes & is.null(attr(ytypes,'CCgroups')) ){
#      groups <- rep(0,S)
#      groups[ytypes == 'CC']  <- 1
#      attr(ytypes,'CCgroups') <- groups
#    }
#    return( list(y = ydata, CCgroups = CCgroups, FCgroups = FCgroups, 
#                 typeNames = ytypes, ydataNames = ydataNames) ) 
#  }
  
  ngroup  <- 0
  ccg     <- CCgroups
  fcg     <- FCgroups
  y       <- numeric(0)

  snames <- colnames(ydata)
  nc     <- ncol(ydata)
  wfact  <- .whichFactor(ydata)
  nfact  <- length(wfact)
  wnot   <- c(1:nc)
  if(nfact > 0)wnot <- wnot[-wfact]
  
  ntypes <- character(0)
  
  if(length(wnot) > 0){
    
    if(is.null(ccg)) ccg <- rep(0,length(wnot)) # if not assigned, assume same 
    if(is.null(fcg)) fcg <- rep(0,length(wnot))
    
    snames <- snames[wnot]
    ntypes <- ytypes[wnot]
    y      <- ydata[,wnot,drop=F]
    
    wcomp <- grep('CC',ytypes[wnot])
    ncomp <- length(wcomp)
    if(ncomp > 0){
      if( max(ccg[wnot[wcomp]]) == 0 )ccg[wnot[wcomp]] <- 1  #assume same group
      goo <- grep('other',snames[wcomp])
      if( length(goo) == 0 )snames[wcomp[ncomp]] <- 
        paste(snames[wcomp[ncomp]],'other',sep='')
    }
    
    wcomp <- grep('FC',ytypes[wnot])
    ncomp <- length(wcomp)
    if(ncomp > 0){
      if( max(fcg[wnot[wcomp]]) == 0)fcg[wnot[wcomp]] <- 1  #assume same group
      goo <- grep('other',snames[wcomp])
      if(length(goo) == 0)snames[wcomp[ncomp]] <- 
        paste(snames[wcomp[ncomp]],'other',sep='')
    }
  }
  
  if(nfact > 0){   # categorical
    
    ngroup <- 0
    ycat   <- cag <- numeric(0)
    if(length(ccg) > 0)cag    <- ccg*0
    
    for(j in 1:nfact){
      
      ngroup <- ngroup + 1
      
      conj   <- contrasts(ydata[,wfact[j]],contrasts=F)
      cj     <- colnames(conj)
      yj     <- conj[ydata[,wfact[j]],]
      colnames(yj) <- paste(colnames(ydata)[wfact[j]],cj,sep='')
      
      w11    <- which(colSums(yj) > 0)  #drop empty levels
      yj     <- yj[,w11]
      cj     <- cj[w11]
      
      goo <- grep('other',colnames(yj))
      if(length(goo) == 0){
        colnames(yj)[ncol(yj)] <- paste(colnames(ydata)[wfact[j]],'other',sep='')
        cj[ncol(yj)] <- colnames(yj)[ncol(yj)]
      }
      
      ycat   <- cbind(ycat, yj)
      cag    <- c(cag,rep(ngroup,length(cj)))
      fcg    <- c(fcg,rep(0,length(cj)))
      ccg    <- c(ccg,rep(0,length(cj)))
      ntypes <- c(ntypes,rep('CAT',length(cj)))
    }
    
    rownames(ycat) <- NULL
    n1 <- ncol(y) + 1
    n2 <- ncol(ycat)
    y <- cbind(y,ycat)
  #  catCols <- c( n1:n2 )
  #  attr(ytypes,'catCols')  <- catCols
    attr(ntypes,'CATgroups')  <- cag
  }
  
  if(max(ccg) > 0)attr(ntypes,'CCgroups') <- ccg
  if(max(fcg) > 0)attr(ntypes,'FCgroups') <- fcg
  
 # y <- as.matrix(y)
 # colnames(y) <- snames

  list(y = as.matrix(y), CCgroups = ccg, FCgroups = fcg, 
       CATgroups = attr(ntypes,'CATgroups'), typeNames = ntypes,
       ydataNames = ydataNames)
}

.setUpSim <- function(n, S, Q, x, typeNames){
  
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  
  notOther <- c(1:S)
  snames   <- character(0)
  tnames   <- character(0)
  sN       <- S
  catCols  <- NULL
  
  ngroup <- fgroup <- cgroup <- 1
  GROUPS <- F
  CCgroups <- FCgroups <- CATgroups <- numeric(0)
  s      <- 0
  
  wcc <- which(!typeNames %in% c('CC','FC','CAT'))
  ncc <- length(wcc)
  if(ncc > 0){
    snames   <- paste('S',c(1:ncc),sep='')
    CCgroups <- FCgroups <- CATgroups <- rep(0,ncc)
    tnames   <- typeNames[wcc]
    s <- ncc
  }
  wcc <- which(typeNames == 'CC')
  ncc <- length(wcc)
  if(ncc > 0){
    ss <- c( (s+1):(s+ncc))
    CCgroups <- c(CCgroups,rep(1,ncc))
    FCgroups <- c(FCgroups,rep(0,ncc))
    tnames   <- c(tnames,rep('CC',ncc))
    snn      <- paste('S',ss,sep='')
    snn[ncc] <- paste(snn[ncc],'other',sep='')
    snames   <- c(snames, snn)
    ngroup <- 1
    s <- max(ss)
  }
  wcc      <- which(typeNames == 'FC')
  ncc <- length(wcc)
  if(ncc > 0){
    ss <- c( (s+1):(s+ncc))
    FCgroups <- c(FCgroups,rep(1,ncc))
    CCgroups <- c(CCgroups,rep(0,ncc))
    tnames   <- c(tnames,rep('FC',ncc))
    snn      <- paste('S',ss,sep='')
    snn[ncc] <- paste(snn[ncc],'other',sep='')
    snames   <- c(snames, snn)
    fgroup <- 1
    s <- max(ss)
  }
  
  CATgroups <- CCgroups*0
  
  if( 'CAT' %in% typeNames ){
    
    wk     <- which(typeNames == 'CAT')
    ncomp  <- length(wk)
    ncat   <- sample(3:4,ncomp,replace=T)
    nall   <- sum(ncat)
    ntot   <- s + nall
    CATgroups <- rep(0,s)
    
    js <- s
    for(j in 1:ncomp){
      js     <- js + 1
      sseq   <- (s+1):(s + ncat[j])
      cj <- paste('S',js,letters[1:ncat[j]],sep='')
      cj[ncat[j]] <- paste('S',js,'other',sep='')
      snames    <- c(snames,cj)
      CATgroups <- c(CATgroups,rep(j,ncat[j]))
      tnames    <- c(tnames,rep('CAT',ncat[j]))
      s <- max(sseq)
    }
    CCgroups <- c(CCgroups,rep(0,sum(ncat)))
    FCgroups <- c(FCgroups,rep(0,sum(ncat)))     
    catCols  <- which(CATgroups > 0)
    cgroup   <- ncomp
  }
  sN     <- length(tnames)
  oo     <- grep('other',snames)
  notOther <- c(1:sN)[-oo]
  
  tmp <- .gjamGetTypes(tnames)
  typeCols <- tmp$typeCols
  typeFull <- tmp$typeFull
  typeCode <- tmp$TYPES[typeCols]
  allTypes <- sort(unique(typeCols))
  typeNames <- tmp$typeNames
  
  if(is.null(x)){
    x <- matrix( rnorm(n*Q,.1), n, Q)  
    x[,1] <- 1
  }
  
  beta <- matrix(0, Q, sN)
  ss   <- diag(.01,sN)    
  
  colnames(beta) <- colnames(ss) <- rownames(ss) <- snames
  wkeep <- numeric(0)
  cnames <- tnames <- character(0)
  
  for(k in allTypes){
    
    wk <- which(typeCols == k)
    nk <- length(wk)
    
    if( typeFull[wk[1]] == 'presenceAbsence' ){
      diag(ss)[wk] <- 1
      beta[,wk]    <- runif(Q*nk,-1.5,1.5)
      wkeep <- c(wkeep,wk)
      tnames <- c(tnames,typeNames[wk])
      cnames <- c(cnames,colnames(beta)[wk])
    }
    if(typeFull[wk[1]] %in% c('continuous','contAbun')){
      diag(ss)[wk] <- .4
      beta[,wk]    <- runif(Q*nk,-.5,2)
      wkeep <- c(wkeep,wk)
      tnames <- c(tnames,typeNames[wk])
      cnames <- c(cnames,colnames(beta)[wk])
    }
    if(typeFull[wk[1]] == 'discAbun'){
      diag(ss)[wk] <- 1
      beta[,wk]    <- runif(Q*nk,-.1,2)
      wkeep <- c(wkeep,wk)
      tnames <- c(tnames,typeNames[wk])
      cnames <- c(cnames,colnames(beta)[wk])
    }
    if(typeFull[wk[1]] == 'ordinal'){
      diag(ss)[wk] <- 1
      beta[,wk]    <- runif(Q*nk,-.4,2)
      wkeep <- c(wkeep,wk)
      tnames <- c(tnames,typeNames[wk])
      cnames <- c(cnames,colnames(beta)[wk])
    }
    
    if( typeFull[wk[1]] %in% c('fracComp','countComp','categorical') ){
      
      if(length(wk) < 2)stop('composition data must have at least 2 columns')
      
      ntt <- cgroup
      if( typeFull[wk[1]] == 'fracComp' ){
        ntt <- fgroup
        attr(tnames,'FCgroups') <- FCgroups
      }
      if( typeFull[wk[1]] == 'countComp' ){
        ntt <- ngroup
        attr(tnames,'CCgroups') <- CCgroups
      }
      if( typeFull[wk[1]] == 'categorical' ){
        attr(tnames,'CATgroups') <- CATgroups
      }
      
      for(i in 1:ntt){
        
        if(ntt == 1){   
          wki <- wk
        } else {
          if( typeFull[wk[1]] == 'countComp' )wki <- 
              which(typeCols == k & CCgroups == i)
          if( typeFull[wk[1]] == 'fracComp' )wki <- 
              which(typeCols == k & FCgroups == i)
          if( typeFull[wk[1]] == 'categorical' )wki <- 
              which(typeCols == k & CATgroups == i)
        }
        
        nki    <- length(wki)
        
        if( typeFull[wk[1]] == 'categorical' ){
          
          bb <- matrix( rnorm(Q*nki,0,.5), Q,nki)
          bb[1,] <- bb[1,]*0
          
          for(kk in 1:5){
            
            mu   <- x%*%bb
            w    <- mu
            cols <- apply(w,1,which.max)
            mindex <- cbind( c(1:n),cols )
            
            wmax <- w[mindex]
            ww   <- which(wmax < 0)
            nw   <- length(ww)
            if(nw > 0) w[mindex[ww,]] <- .tnorm(nw,0,10,mu[mindex[ww,]],1)
            
            bb <- solve(crossprod(x))%*%crossprod(x,w)
          }
          
          keep <- as.numeric( names(table(cols))  )
          wkeep <- c(wkeep,wki[keep])
          tnames <- c(tnames,rep('CAT',length(keep)))
          
          bbb  <- colnames(beta)[wki[keep]]
          if(length(keep) < nki){
            bbb  <- substr(bbb,1,2)
            labs <- c(letters[1:(length(bbb) - 1)],'other')
            bbb  <- paste(bbb,labs,sep='')
          }
          cnames <- c(cnames,bbb)
          beta[,wki] <- bb
          diag(ss)[wk] <- 1
          
        } else {
          
          bb <- matrix( rnorm(Q*nki,0,1/nki), Q, nki)
          bb[1,] <- bb[1,]*0
            
            w  <- x%*%bb
            
            for(m in 1:3){
              w1 <- w
              w1[w < 0] <- 0
              w2 <- sweep(w1,1,rowSums(w1),'/')
              w[w >= 0] <- w2[w >= 0]
              bb <- solve(crossprod(x))%*%crossprod(x,w)
              w  <- x%*%bb
            }

          wkeep <- c(wkeep,wki)
          tnames <- c(tnames,typeNames[wki])
          cnames <- c(cnames,colnames(beta)[wki])
          diag(ss)[wk] <- .1/nk^2.5
          beta[,wki] <- bb
        }
        
      }
    }
    
  }
  
  S <- length(wkeep)
  beta      <- beta[,wkeep]
  sigma     <- ss[wkeep,wkeep]
  colnames(beta) <- colnames(sigma) <- rownames(sigma) <- cnames
  CCgroups  <-  CCgroups[wkeep] 
  FCgroups  <-  FCgroups[wkeep]
  CATgroups <-  CATgroups[wkeep]
  snames    <- cnames
  other     <- numeric(0)
  notOther  <- c(1:S)
  other     <- grep('other',snames)
  if(length(other) > 0)notOther <- notOther[-other]
  
  list(beta = beta, x = x, sigma = sigma, CCgroups = CCgroups, 
       FCgroups = FCgroups, CATgroups = CATgroups, typeNames = tnames,
       other = other, notOther = notOther, snames = snames)
}

.between <- function(x,lo,hi,ILO = T, IHI = T, OUT=F){
  
  if(length(x) == 0) return( numeric(0) )
  
  if(OUT)return( which(x < lo | x > hi) )
  if(!ILO & !IHI ) return( which(x > lo & x < hi) )
  if(!ILO &  IHI ) return( which(x > lo & x <= hi) )
  if( ILO & !IHI ) return( which(x >= lo & x < hi) )
  if( ILO &  IHI ) return( which(x >= lo & x <= hi) )
        
}        

.simData <- function( n, S, Q, x, typeNames, nmiss, effort ){
  
#  pg <- .95
  
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  
  typeNotCat <- typeNames
  
  cgrep <- grep('CAT',typeNames)
  if(length(cgrep) > 0){
    ycat <- vector( mode = 'list', length=length(cgrep) )
    names(ycat) <- paste('CAT',1:length(cgrep),sep='_')
  }
    
  cuts <- numeric(0)
  
  tmp    <- .setUpSim(n, S, Q, x, typeNames)
  beta   <- tmp$beta
  x      <- tmp$x
  sig    <- tmp$sigma
  snames <- colnames(beta)
  typeNames <- tmp$typeNames
  other     <- tmp$other
  notOther  <- tmp$notOther
  CCgroups  <- tmp$CCgroups
  FCgroups  <- tmp$FCgroups
  CATgroups <- tmp$CATgroups
      
  tmp <- .gjamGetTypes(typeNames)
  typeCols  <- tmp$typeCols
  typeFull  <- tmp$typeFull
  typeCode  <- tmp$TYPES[typeCols]
  allTypes  <- sort(unique(typeCols))
  
  S      <- length(typeNames)
  xnames <- paste('x',1:Q,sep='')
  
  SS    <- matrix(1,S,S)
  SS[lower.tri(SS)] <- runif(S*(S - 1)/2,-.98,.98)
  SS[upper.tri(SS)] <- SS[lower.tri(SS)]
  
  SS    <- cor( .rMVN(S+5,0,SS) )
  SS    <- .cor2Cov(diag(sig),SS)
  
  sigma <- .rwish(S+2,SS)/(S + 2)
  
  corCols <- which(typeNames %in% c('PA','OC','CAT'))
  if(length(corCols) > 0){
    corSpec <- .cov2Cor(sigma)
    sigma[corCols,corCols] <- corSpec[corCols,corCols]
  }
  
  beta[,other] <- 0
  mu <- w <- matrix(0,n,S)
  
  mu[,notOther] <- x%*%beta[,notOther]
  w[,notOther]  <- mu[,notOther] + .rMVN(n,0,sigma[notOther,notOther]) 
  colnames(w) <- snames
    
  y  <- w
  z  <- w*0
  z[w <= 0]   <- 1
  z[w > 0]    <- 2
  
  for(k in allTypes){
    
    wk <- which(typeCols == k)
    nk <- length(wk) 
    
    if( typeFull[wk[1]] %in% c('fracComp','countComp','categorical') ){
      
      if( typeFull[wk[1]] == 'fracComp' )
        groups <- attr(typeNames,'FCgroups') <- FCgroups
      if( typeFull[wk[1]] == 'countComp' )
        groups <- attr(typeNames,'CCgroups') <- CCgroups
      if( typeFull[wk[1]] == 'categorical' )
        groups <- attr(typeNames,'CATgroups') <- CATgroups
      
      ntt <- max(c(1,groups))
      
      for(i in 1:ntt){
        
        if(ntt == 1){
          wki <- wk
        } else {
          wki <- which(typeCols == k & groups == i)
        }
        nki <- length(wki)
        
        if( typeFull[wk[1]] == 'categorical' ){
          
          wko  <- wki[1:(nki-1)]                  
          wcol <- apply(w[,wko],1,which.max)
          w0   <- which( w[,wko][ cbind( c(1:n),wcol ) ] < 0 )
          if(length(w0) > 0)wcol[w0] <- nki
          
          wtab <- tabulate(wcol)
          if(length(wtab) < nki){
            ww <- rep(0,nki)
            ww[1:length(wtab)] <- wtab
            wtab <- ww
          }
          
          if(min(wtab) < 5){
            wlo <- which(wtab < 5)
            for(s in 1:length(wlo)){
              wro <- sample(n,5)
              wcol[wro] <- wlo[s]
              tmp <- w[wro,wki]
              if(wlo[s] == nki){
                tmp[tmp > -.01] <- -.01   # all values neg
                tmp[,nki] <- .1
              } else {
                mm <- pmax(0,apply(tmp,1,max))
                tmp[,wlo[s]] <- mm + .1
              }
              w[wro,wki] <- tmp
            }
          }
            
          mindex <- cbind(1:n,wcol)
          
          vv <- colnames(w)[wki[wcol]]
          mm <- nchar(vv)
          vv <- substr(vv,3,mm)                  
          
          ycat[[i]] <- vv
          
          yk   <- w[,wki]*0
          yk[ mindex ] <- 1
          y[,wki] <- yk
          z[,wki] <- yk + 1
          
        } else {
          
          noto <- c(1:nki)[-grep('other',snames[wki])]
          
          ww     <- w[,wki]
          
          for(j in 1:5){
            
            w0     <- which(ww < 0)
            ww[w0] <- 0  
            
            yk  <- .gjamCompW2Y(ww,notOther=noto)$ww
            
            yplus <- which(yk > 0)
            yminu <- which(yk < 0)
            
            ww[yplus] <- yk[yplus]
            
            bb <- solve(crossprod(x))%*%crossprod(x,ww)
            mu <- x%*%bb
            ww <- mu + .rMVN(n,0,sigma)[,wki]
          }
          zk     <- ww*0 + 1
          zk[w0] <- 0
          w[,wki] <- ww  
          beta[,wki] <- bb
          
          if(typeFull[wk[1]] == 'fracComp'){
            y[,wki] <- yk
            z[,wki] <- zk
          }
          if( typeFull[wk[1]] == 'countComp' ){
            
            mm <- S*20
            a  <- 4
            b  <- mm/a
            ee <- rpois(n,rgamma(n,shape=a,scale=b))
            yy <- sweep(yk,1,ee,'*')
            
            ww <- ceiling(yy)
            ww[ww < 0] <- 0
            y[,wki] <- ww
            z[,wki] <- ww + 1
          }
        }
      }
    }
    
    if( typeFull[wk[1]] != 'continuous' ) y[,wk][y[,wk] < 0] <- 0  # not cens
    
    if( typeFull[wk[1]] == 'presenceAbsence' )y[,wk] <- z[,wk] - 1       
    
    if( typeFull[wk[1]] == 'discAbun' ){
      
      if(!is.null(effort)){
        we     <- wk[wk %in% effort$columns]
        y[,we] <- round( w[,we]*effort$values,0 )
      } else {
        w0 <- round(w[,wk,drop=F],0)
        y[,wk] <- w0
      }
      y[,wk][y[,wk] < 0] <- 0
      z[,wk]        <- y[,wk] + 1
    }
    
    if( typeFull[wk[1]] == 'ordinal' ){
      
      yy   <- w[,wk,drop=F]
      ncut <- 8
      maxw <- floor(max(yy))
      
      cuts  <- t( matrix( c(-Inf, seq(0,(maxw-1),length=(ncut-2)) ,Inf),
                          ncut,nk) )
      rownames(cuts) <- snames[wk]
      
      for(j in 1:nk){
        z[,wk[j]]   <- findInterval(yy[,j],cuts[j,])
      }
      
      y[,wk] <- z[,wk] - 1
  #    cuts   <- .gjamTheta2cuts(cuts,sigma[wk,wk])   
    }
    
  }
  
  #####################################
  
  noMore <- F
  if( 'categorical' %in% typeFull & noMore){
    
  #  wss  <- w
  #  w[,notOther] <- .sqrtRootMatrix(w[,notOther],sg[notOther,notOther],DIVIDE=T)
  #  css  <- .cov2Cor(sg[notOther,notOther])
  #  mus <- x%*%beta
    wss  <- w*0
    wss[,notOther] <- .sqrtRootMatrix(w[,notOther],sigma[notOther,notOther],
                                      DIVIDE=T)
    css  <- .cov2Cor(sigma[notOther,notOther])
    alpha <- .sqrtRootMatrix(beta,sigma, DIVIDE=T)
    muss <- x%*%alpha
    
    wk <- which(typeNames == 'CAT')
 #   nk <- length(wk)
    wo <- which(wk %in% notOther)
 #   wu <- which(typeCols[notOther] == k)
 #   wp <- w[, wk, drop=F]*0
    
    plo  <- w*0 - 500
    phi  <- w*0 + 500
    
    phi[y == 0] <- 0
    plo[y == 1] <- w[y == 1]
    IXX <- solve(crossprod(x))
    
    for(k in 1:25){
     
      tmp <- .gjamWcatLoop2(y, ws = wss, mus = muss, sgs = css, 
                            notOther, plo, phi, groups = CATgroups)
      #     tmp <- .gjamWcatLoop2(y, ws = w, mus = muw, sgs = sg[notOther,notOther], 
      #                           notOther, plo, phi, groups = CATgroups)
      wss[,wk] <- tmp$w[,wk]
      plo     <- tmp$plo
      phi     <- tmp$phi
      beta[,wo] <- IXX%*%crossprod(x,wss[,wo])
      muss[,wo] <- x%*%beta[,wo]
    }
    w[,wo] <- wss[,wo]
  }
      
      
  beta <- solve(crossprod(x))%*%crossprod(x,w)
  sigma[notOther,notOther] <- var(w[,notOther] - x%*%beta[,notOther]) ### NO
  sigma[other,] <- sigma[,other] <- 0
  diag(sigma)[other] <- diag(sig)[other]
  
  ydata <- data.frame(y)
  typeFrame <- typeNames
  
  if('CAT' %in% typeNames){
    wcat <- grep('CAT',typeNames)
    wnot <- c(1:S)[-wcat] 
    nss  <- length(wnot) + 1
    ncc  <- length(wnot) + length(ycat)
    names(ycat) <- paste('S',nss:ncc,sep='')
    ydata <- as.data.frame(ycat)
    if(length(wnot) > 0)ydata <- cbind(y[,wnot,drop=F],ydata)
    typeFrame <- c(typeNames[wnot], rep('CAT',length(ycat)))
  }
  
  if(nmiss > 0){
    x[ sample(length(x),nmiss) ] <- NA
    x[,1] <- 1
    wmiss <- which(is.na(x),arr.ind=T)
    nmiss <- nrow(wmiss)
  }
  
  xnames[1]      <- 'intercept'
  colnames(y)    <- snames
  colnames(beta) <- rownames(sigma) <- colnames(sigma) <- snames
  colnames(x)    <- rownames(beta) <- xnames
  
  form <- as.formula( paste('~ ',paste(colnames(x)[-1],collapse='+' )) )

  list(formula = form, xdata = data.frame(x), ydata = ydata,
       y = y, w = w,  typeNames = typeFrame, typeY = typeNames, effort = effort,
       trueValues = list(beta = beta, sigma = sigma, 
                         corSpec = .cov2Cor(sigma), cuts = cuts))
}

.tnorm <- function(n,lo,hi,mu,sig){   

  #normal truncated lo and hi
  
  tiny <- 10e-6

  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))

  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 

  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  z[z == Inf]  <- lo[z == Inf] + tiny
  z[z == -Inf] <- hi[z == -Inf] - tiny
  z
}
.traitLabel <-
function(tname){
  
  tname <- .replaceString(tname,now='soilFactor',new='')
  tname[tname == 'gmPerSeed'] <- 'Seed mass'
  tname[tname == 'gmPerCm']   <- 'Wood dens'
  tname[tname == 'woodSG']    <- 'Wood dens (green)'
  tname[tname == 'maxHt']     <- 'Max ht'
  tname[tname == 'leafN']     <- 'leaf [N]'
  tname[tname == 'leafP']     <- 'leaf [P]'
  tname[tname == "other"]  <- 'Deciduous'
  tname[tname == "broaddeciduous"]  <- 'Deciduous'
  tname[tname == "broadevergreen"]  <- 'BL evergrn'
  tname[tname == "needleevergreen"] <- 'NL evergrn'
  tname[tname == "dioecious"] <- 'Dioecious'
  tname[tname == "u1"] <- 'Slope'
  tname[tname == "u2"] <- 'Aspect 1'
  tname[tname == "u3"] <- 'Aspect 2'
  tname[tname == "ringPorous"] <- 'RP xylem'
  tname[tname == "temp"] <- 'Winter temperature'
  tname[tname == "stdage"] <- 'Stand age'
  for(j in length(tname)){
    tname[j] <- paste(toupper(substring(tname[j], 1, 1)), substring(tname[j], 2),sep = "", collapse = " ")
  }
  tname
}
.updateWishartNoPrior <- function(xx,yy,df,beta=NULL,IXX=NULL,WX=NULL,WIX=NULL,
                                  TRYPRIOR=F){
  
  # more stable without prior
  # TRYPRIOR includes non-informative prior if cholesky fails 
  
  index <- 0
  
  if(is.null(IXX)) IXX <- chol2inv(chol( crossprod(xx) ) )
  if(is.null(WX))  WX  <- crossprod(xx,yy)
  if(is.null(WIX)) WIX <- IXX%*%WX
  
  SSS   <- crossprod(yy) - t(WX)%*%WIX
  testv <- try(chol(SSS),T)
  
  if( inherits(testv,'try-error') ){
    tiny <- 1e-8
    SSS[SSS < tiny] <- tiny
    message('warning: updateWishartNoPrior')
    SSS <- crossprod(yy - xx%*%beta) +  diag(diag(SSS)*.001)#*nrow(SSS)
    SSS <- SSS + diag(diag(SSS)*.1)
    testv <- try(chol(SSS),T)
  }
  
  SI <- chol2inv(testv)
  
  testChol <- try(chol(SI),T)
    
  if( inherits(testChol,'try-error') ){
    message('warning: prior used in updateWishartNoPrior')
    if(TRYPRIOR){
      index  <- 1
      SI     <- SI + diag(diag(SI)*.01)
      df     <- nrow(SI) + nrow(xx)
      testChol <- try(chol(SI),T)
    }
  }

  z     <- matrix(rnorm(df*nrow(SSS)),df,nrow(SSS))%*%testChol
  sinv  <- crossprod(z)
  
  testSolve <- try( solve(sinv),T )
  if( !inherits(testSolve,'try-error') )sigma <- testSolve
  
  if( inherits(testSolve,'try-error') ){
    message('warning: prior used in updateWishartNoPrior')
    if(TRYPRIOR){
      sinv   <- sinv + diag(diag(sinv)*.01)
      df     <- nrow(sinv) + nrow(xx)
      testSolve <- try(chol(sinv),T)
      sigma   <- chol2inv(testSolve)
    }
  }
  list( sigma = sigma, sinv = sinv, indicator = index )
}

.sqrtRootMatrix <- function(xmat,sigma,DIVIDE=F){
  
  # xmat is n by p
  # sigma is p by p
  
  if(DIVIDE){
    if(length(sigma) == 1)return(xmat/sqrt(sigma))
    return( xmat%*%diag(1/sqrt(diag(sigma))) )
  }
  
  if(length(sigma) == 1)return(xmat*sqrt(sigma))
  xmat%*%diag(sqrt(diag(sigma)) )
}

.yaxisHorizLabs <- function( labels, at=c(1:length(labels)), xshift=.05,
                             col = 'black', pos=NULL){
  
  #add horizontal y axis labels to existing plot
  #pos should be either NULL, 2 (left)
  
  text(par('usr')[3] - xshift*par('usr')[4] - par('usr')[3], y=at,
       labels, xpd=T, pos = pos, col=col)
}
###################################
.sample.p <- function(N,avec,bvec,K){
  
  a    <- avec + vapply(1:(N-1),function(k)sum(K==k), 0)
  b    <- bvec + vapply(1:(N-1),function(k)sum(K>k), 0)
  V    <- rbeta((N-1), a, b)
  p    <- vector("numeric",length=N)
  p[1] <- V[1]
  for(l in 2:(N-1))p[l] <- prod(1 - V[1:(l-1)])*V[l]
  p[N] <- prod(1 - V)   
  p
}

.getPars <- function(x, N, r, Y, B, D, Z, sigmaerror, K, pvec,
                     alpha.DP, inSamples,lo=NULL,hi=NULL,...){      
  ###r = 5, alpha.DP=20
  
  X  <- x[inSamples,]
  XX <- crossprod(X)
  
  if(is.null(lo) & is.null(hi)){
    .sample.B <- function(Y,W,A,sigma.sq){
      OmegaB  <- .solveArma(((1/sigma.sq)*XX))
      muB     <- t(OmegaB%*%crossprod((1/sigma.sq)*X,(Y - W%*%t(A))))
      B       <- .rmvnormArma(ncol(Y),rep(0,ncol(OmegaB)),OmegaB) + muB
      return(B)
    }
  } else{
    .sample.B <- function(Y,W,A,sigma.sq){
      OmegaB  <- .solveArma(((1/sigma.sq)*XX))
      muB     <- t(OmegaB%*%crossprod((1/sigma.sq)*X,(Y - W%*%t(A))))
      B <- .tnormMVNmatrixRcpp(avec=muB,muvec=muB,smat=OmegaB,
                                 lo=lo,hi=hi)
      return(B)
    }
  }
  
  nn   <- length(inSamples)
  ntot <- nrow(Y)
  p    <- ncol(x)
  S    <- ncol(Y)
  
  #------ rnd effects for the variance matrix "R", step 2
  covR <- .solveArma( (1/sigmaerror)*crossprod(Z[K,]) + diag(r) ) # Sigma_W
  z1   <- crossprod( Z[K,]/sigmaerror,t(Y - x%*%t(B)) )        
  RR   <- .rmvnormArma(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
  if(nn < ntot)RR[-inSamples,] <- .rmvnormArma(ntot-nn,mu=rep(0,r),
                                               sigma=diag(r))
  
  #------ potential vector values Z
  avec <- 1/rgamma(r, shape = (2 + r )/2, 
                    rate = ((1/1000000) + 2*diag(.solveArma(D)) ) )  
  
  D    <- .riwish(df = (2 + r + N - 1), S = (crossprod(Z) + 2*2*diag(1/avec)))
  Z    <- .bodyfnZarma(kk=K, Yk=Y[inSamples,], Xk=x[inSamples,], Dk=D,Bk=B, 
                       Wk=RR[inSamples,], sigmasqk=sigmaerror, Nz=N)
  
  #------ sigma error
  RndEff <- RR%*%t(Z[K,])
  res    <- sum((Y[inSamples,] - x[inSamples,]%*%t(B) - 
                   RndEff[inSamples,] )^2)
  sigmaerror <- 1/rgamma(1,shape=(S*nn + 1)/2,rate=res/2)  # no prior?
  
  #------ betas
  A <- Z[K,]
  B <- .sample.B(Y = Y[inSamples,], W = RR[inSamples,], A = A, 
                 sigma.sq = sigmaerror)
  
  #------ labels "K"
  
 # ytmp  <- .sqrtRootMatrix(Y[inSamples,],sg, DIVIDE=T)    # new Yk
 # btmp  <- .sqrtRootMatrix(t(B),sg)             # new B
                          
  pmat <- .obtainpmatKcpp(pveck = pvec,Yk = Y[inSamples,], Zk = Z,
                          Xk = x[inSamples,], Bk = B, Wk = RR[inSamples,],
                          sigmasqk = sigmaerror)
  K <- unlist( apply(pmat, 1, function(x)sample(1:N,size=1,prob=x)) )
  
  #------ cluster probs "pvec"
  pvec <- .sample.p(N = N, avec = rep(alpha.DP/N,(N-1)),
                    bvec = ((N-1):1)*alpha.DP/N, K = K)  
  
  list(A = A, D = D, Z = Z, B = B, K = K, pvec = pvec, 
       sigmaerror = sigmaerror, RndEff = RndEff)
} 

.wWrapper <- function(REDUCT, n, S, effMat, corCols, typeNames, 
                      typeFull, typeCols, 
                      allTypes, holdoutN, holdoutIndex, censor, 
                      censorCA, censorDA, notOther, sampleW, byRow, byCol,
                      indexW, ploHold, phiHold, sampleWhold)
  if(REDUCT){
    
    function(x, w, y, muw, sg, alpha, cutg, plo, phi, rndEff, sigmaerror, wHold){
      
      SC   <- ncol(y)
      scol <- c(1:S)
      muf <- muw + rndEff
      
      if(holdoutN > 0)wHold <- w[drop=F,holdoutIndex,]
      
      if(length(corCols) > 0){
        scol <- scol[-corCols]
        SC   <- length(scol)
      }
      
      yPredict  <-  w*0
      
      mef <- as.vector(muf[,scol])
      w[,scol]  <- matrix( .tnorm(n*SC, as.vector(plo[,scol]), 
                                  as.vector(phi[,scol]),mef,
                                  sqrt(sigmaerror)),n,SC)        # cov scale
      yPredict[,scol] <- matrix( rnorm(n*SC,mef,sqrt(sigmaerror)),n,SC)
      
      if(holdoutN > 0){
        wHold[,scol] <- matrix( .tnorm(holdoutN*SC, as.vector(ploHold[,scol]), 
                                           as.vector(phiHold[,scol]),
                                           as.vector( muf[drop=F,holdoutIndex,scol] ),
                                           sqrt(sigmaerror)),holdoutN,SC) 
      }
      
      sigvec <- rep(sigmaerror,S)
      
      if(length(corCols) > 0){   # corr scale (spp have different sigmaerror)
        SC     <- length(corCols)
        mef    <- .sqrtRootMatrix(muw[,corCols] + rndEff[,corCols],
                                  sg[corCols,corCols], DIVIDE=T)
        
        mef    <- as.vector(t(mef))
        w[,corCols] <- matrix( .tnorm(n*SC, as.vector(t(plo[,corCols])), 
                                       as.vector(t(phi[,corCols])), mef,1),
                                       n,SC, byrow=T) # cor scale
        if(holdoutN > 0){
          wHold[,corCols] <- matrix( .tnorm(holdoutN*SC, as.vector(ploHold[,corCols]), 
                                             as.vector(phiHold[,corCols]),mef,
                                             sqrt(sigmaerror)),holdoutN,SC) 
        }
        yPredict[,corCols] <- matrix( rnorm(n*SC, mef,1),n,SC, byrow=T)
      }
      
      
      w[sampleW == 0] <- y[sampleW == 0]
      if(holdoutN > 0){
        wHold[sampleWhold == 0] <- y[holdoutIndex,][sampleWhold == 0]
      }
        
      FCgroups  <- attr(typeNames,'FCgroups')
      CCgroups  <- attr(typeNames,'CCgroups')
      CATgroups <- attr(typeNames,'CATgroups')
      
      for(k in allTypes){
        
        wk <- which(typeCols == k)
        wo <- which(wk %in% notOther)
        nk <- length(wk)
        wu <- which(typeCols[notOther] == k)
        wp <- w[, wk, drop=F]
        yp <- yPredict[, wk, drop=F]
        
        groups <- NULL
        if(typeFull[wk[1]] == 'countComp')  groups <- CCgroups[wk]
        if(typeFull[wk[1]] == 'fracComp')   groups <- FCgroups[wk]
        
        #################
        
        if( typeFull[wk[1]] == 'categorical' ){
          groups <- CATgroups[wk]
          tmp <- .gjamWcatLoop2(y, ws = wp, mus = muf, sgs = sigvec, 
                                notOther = notOther, plo, phi, 
                                groups = CATgroups, REDUCT=T)
          wp[,wo] <- tmp$w[,wo]
          plo     <- tmp$plo
          phi     <- tmp$phi
          
          if(holdoutN > 0){
            wHold[,wo] <- .gjamWcatLoop2(y, ws = w[, wk, drop=F], 
                                                  mus = muf, sgs = sigvec, 
                                                  notOther = notOther, ploHold, phiHold, 
                                                  groups = CATgroups, REDUCT=T) 
          }
        }
        
        tmp <- .gjamWLoopTypes(wo, typeFull[wk[1]], y[,wk,drop=F], wp, yp, cutg, 
                               censor, censorCA, censorDA, effMat[,wk,drop=F], 
                               groups, k, typeCols, notOther, wk = wk)
        w[,wk]        <- tmp[[1]]
        yPredict[,wk] <- tmp[[2]]
        
        if(holdoutN > 0){
          wp  <- wHold[,wk,drop=F]
          yp  <- yPredict[holdoutIndex, wk, drop=F]
   #       eff <- list(columns = effort$columns, 
   #                   values = effort$values[drop=F,holdoutIndex, ])
          tmp <- .gjamWLoopTypes(wo, typeFull[wk[1]], y[holdoutIndex,wk,drop=F], wp, 
                                 yp, cutg, censor, censorCA, censorDA, 
                                 effMat[holdoutIndex, wk, drop=F], 
                                 groups, k, typeCols, notOther, wk = wk)
          wHold[,wk] <- tmp[[1]]
          yPredict[holdoutIndex,wk] <- tmp[[2]]
        }
        yPredict[,wk] <- .censorValues(censor,y,yPredict)[,wk]
      }
      
      
      w[sampleW == 0] <- y[sampleW == 0]
      if(holdoutN > 0){
        wHold[sampleWhold == 0] <- y[holdoutIndex,][sampleWhold == 0]
      }
      
      list(w = w, wHold = wHold, yp = yPredict, plo = plo, phi = phi )
    }
    
  } else {
    
    function(x, w, y, muw, sg, alpha, cutg, plo, phi, rndEff  = NULL, 
             sigmaerror = NULL, wHold){
      
      w[sampleW == 0] <- y[sampleW == 0]
      if(holdoutN > 0){
        wHold[sampleWhold == 0] <- y[holdoutIndex,][sampleWhold == 0]
      }
      
      yPredict <- w*0
      yPredict[,notOther] <- .rMVN(n,muw[,notOther],sg[notOther,notOther])
      ypred <- yPredict
      
      if( length(corCols) > 0 ){    #expanded w on this scale
        wss  <- w*0
        css  <- .cov2Cor(sg[notOther,notOther])
        muss <- x%*%alpha
        ypred[,notOther]   <- .rMVN(n,muss[,notOther],css)
        yPredict[,corCols] <- ypred[,corCols]
      } 
      FCgroups  <- attr(typeNames,'FCgroups')
      CCgroups  <- attr(typeNames,'CCgroups')
      CATgroups <- attr(typeNames,'CATgroups')
      
      for(k in allTypes){
        
        wk <- which(typeCols == k)
        nk <- length(wk)
        wo <- which(wk %in% notOther)
        wu <- which(typeCols[notOther] == k)
        wp <- w[, wk, drop=F]*0
        yp <- yPredict[, wk, drop=F]
        
        if( typeFull[wk[1]] %in% c('presenceAbsence','ordinal') ) {
          wss[,notOther] <- .sqrtRootMatrix(w[,notOther],sg[notOther,notOther],
                                            DIVIDE=T)               ###NOT
          wp[,wo] <- .gjamWLoop(ws = wss[,notOther], mus = muss[,notOther], 
                                sgs = css, wkk = wu, 
                                lo = plo[,notOther], hi = phi[,notOther],
                                sampW = sampleW[,notOther], indexW)[,wu]
          
          if(holdoutN > 0){
            wHold[,wo] <- .gjamWLoop(ws = wss[drop=F,holdoutIndex,notOther], 
                                              mus = muss[drop=F,holdoutIndex,notOther], 
                                              sgs = css, wkk = wu, 
                                              lo = ploHold[drop=F,,notOther], 
                                              hi = phiHold[drop=F,,notOther],
                                              sampW = sampleWhold[,notOther], 
                                              indexW)[,wu] 
          }
        }
        
        if( !typeFull[wk[1]] %in% c('presenceAbsence','ordinal','categorical') ){
          
          wp[,wo] <- .gjamWLoop(w[,notOther], mus = muw[,notOther], 
                                sgs = sg[notOther,notOther], wkk = wu,
                                plo[,notOther],phi[,notOther],
                                sampleW[,notOther], indexW, byCol, byRow)[,wu]
          if(holdoutN > 0){
            wHold[,wo] <- .gjamWLoop(w[drop=F,holdoutIndex,notOther], 
                                              mus = muw[drop=F,holdoutIndex,notOther], 
                                              sgs = sg[notOther,notOther], wkk = wu,
                                              ploHold[drop=F,,notOther],
                                              phiHold[drop=F,,notOther],
                                              sampleWhold[,notOther], indexW=wo, 
                                              byCol, byRow)[,wu] 
          }
        }
        
        if( typeFull[wk[1]] == 'categorical' ){
          wss[,notOther] <- .sqrtRootMatrix(w[,notOther],sg[notOther,notOther],
                                            DIVIDE=T)
          yy  <- y
          if(holdoutN > 0)yy[holdoutIndex,] <- yp[holdoutIndex,]
          tmp <- .gjamWcatLoop2(yy, ws = wss, mus = muss, sgs = css, 
                                notOther, plo, phi, groups = CATgroups)
          wp      <- tmp$w[,wk]
          plo     <- tmp$plo
          phi     <- tmp$phi
          
          if(holdoutN > 0){
            wHold[,wk] <- .gjamWcatLoop2(yp[drop=F,holdoutIndex,],
                                                  wss[drop=F,holdoutIndex,], 
                                                  muss[drop=F,holdoutIndex,], sgs = css, 
                                                  notOther, ploHold, phiHold, 
                                                  groups = CATgroups)$w[,wk]
          }
        }
        
    #    if(holdoutN > 0)wp[holdoutIndex,wo] <- yp[holdoutIndex,wo]
        
        groups <- NULL
        if(typeFull[wk[1]] == 'countComp')  groups <- CCgroups[wk]
        if(typeFull[wk[1]] == 'fracComp')   groups <- FCgroups[wk]
        if(typeFull[wk[1]] == 'categorical')groups <- CATgroups[wk]
        
        tmp <- .gjamWLoopTypes(wo, type = typeFull[wk[1]], yy = y[,wk,drop=F], 
                               wp, yp, cutg, 
                               censor, censorCA, censorDA, effMat[,wk,drop=F], groups, 
                               k, typeCols, notOther, wk = wk )
        w[,wk]        <- tmp[[1]]
        
        if(holdoutN > 0){
          w[holdoutIndex,wk] <- .gjamWLoopTypes(wo, type = typeFull[wk[1]], 
                                                 yy = yp[holdoutIndex,,drop=F], 
                                                 wp[drop=F,holdoutIndex,],
                                                 yp[drop=F,holdoutIndex,], cutg, 
                                                 censor, censorCA, censorDA, 
                                                 effMat[drop=F,holdoutIndex,], groups, 
                                                 k, typeCols, notOther, wk = wk )[[1]]
        }
        
        yPredict[,wk] <- tmp[[2]]
        yPredict[,wk] <- .censorValues(censor, y, yPredict)[,wk]
      }
      
      w[sampleW == 0] <- y[sampleW == 0]
      if(holdoutN > 0){
        wHold[sampleWhold == 0] <- y[holdoutIndex,][sampleWhold == 0]
      }
      
      list(w = w, wHold = wHold, yp = yPredict, plo = plo, phi = phi )
    }
  }

.paramWrapper <- function(REDUCT,inSamples,SS,loB,hiB,updateBeta){   
  
  if(REDUCT){    
    
    function(x,beta,Y,otherpar){
      
      N  <- otherpar$N
      r  <- otherpar$r
      D  <- otherpar$D
      Z  <- otherpar$Z
      sigmaerror <- otherpar$sigmaerror
      K          <- otherpar$K
      pvec       <- otherpar$pvec
      alpha.DP   <- otherpar$alpha.DP
      lo <- hi <- NULL
      if(!is.null(loB))lo <- t(loB)
      if(!is.null(hiB))hi <- t(hiB)
      tmp        <- .getPars(x = x, N = N, r = r, Y = Y, B = t(beta), 
                             D = D, Z = Z,
                            sigmaerror = sigmaerror,
                            K = K, pvec = pvec, alpha.DP = alpha.DP,
                            inSamples = inSamples, SELECT = F, 
                            SelPars = otherpar$SelPars, lo=lo, hi=hi)
      
      sg <- with(tmp, .expandSigma(sigma = tmp$sigmaerror, SS, Z = tmp$Z, 
                                   K = tmp$K, REDUCT=T))
      
      otherpar=list(A = tmp$A, N = N, r = r, D = tmp$D, Z = tmp$Z, 
                    sigmaerror = tmp$sigmaerror,
                    pvec = tmp$pvec, K = tmp$K, alpha.DP = alpha.DP,
                    SelPars = tmp$SelPars)
      
      return(list(sg = sg, bg = t(tmp$B), rndEff = tmp$RndEff, otherpar = otherpar))
    }
    
  } else {
    
    function(x,beta,Y,otherpar){
      
      sigmaDf  <- otherpar$sigmaDf
      XX  <- crossprod(x[inSamples,])
      IXX <- solve(XX)
      WX  <- crossprod(x[inSamples,],Y[inSamples,])
      WIX <- IXX%*%WX
      
      sg <- .updateWishartNoPrior( x[inSamples,],Y[inSamples,],sigmaDf,
                                   beta=beta,IXX=IXX,WX=WX,WIX=WIX,
                                   TRYPRIOR=T)$sigma
      tmp <- updateBeta(WIX=WIX, IX=IXX, sg=sg, w=Y[inSamples,],
                        y0=Y, inSamples=inSamples,
                        alpha=beta, loBeta=loB, hiBeta=hiB)
      
      bg <- alpha <- tmp$bg
      
      otherpar=list(Z=NA,K=NA,sigmaDf=sigmaDf)
      
      return(list(sg = sg, bg = bg, otherpar = otherpar))
    }
  }
}

.rwish <- function(df,SS){
  z  <- matrix(rnorm(df*nrow(SS)),df,nrow(SS))%*%chol(SS)
  crossprod(z)
}

.riwish <- function(df,S){
  solve(.rwish(df,solve(S)))
}

.expandSigmaChains <- function(snames, sgibbs, otherpar, 
                               simIndex = sample(nrow(sgibbs),50,replace=T), 
                               sigErrGibbs, kgibbs=NULL, 
                               REDUCT=F, CHAINSONLY=F){
  tiny <- 1e-8
  
  S <- otherpar$S
  K <- otherpar$K
  N <- otherpar$N
  r <- otherpar$r
  if(length(simIndex) > 1000)simIndex <- sample(simIndex,1000)
  ns     <- length(simIndex)
  xnames <- otherpar$xnames
  
  if(CHAINSONLY & !REDUCT){  #only return expanded sgibbs
    
    imat   <- matrix(1:(S*S),S,S)
    jmat   <- matrix(1:(S*S),S,S,byrow=T)
    tmp    <- matrix(NA,nrow(sgibbs),S*S)
    sindex <- imat[lower.tri(imat,diag=T)]
    tmp[,sindex] <- sgibbs
    sindex <- jmat[lower.tri(imat,diag=T)]
    tmp[,sindex] <- sgibbs
    
    sMu <- matrix( colMeans(tmp),S,S)
    sSe <- matrix( apply(tmp,2,sd),S,S)
    
    chainList <- list(cchain = NULL, schain = tmp, kchain = NULL)
    
    return( list(chainList = chainList, rMu = NULL, rSe = NULL, 
                 sMu = sMu, sSe = sSe) )
  }
  
  # summarize chains
  
  other    <- grep('other',snames)
  notOther <- c(1:S)
  if(length(other) > 0)notOther <- notOther[-other]
  
  Kindex <- which(lower.tri( diag(S),diag=T ) )
  kchain <- NULL
  
  schain <- cchain <- matrix(0,ns,length(Kindex))
  if(REDUCT)kchain <- matrix(0,ns,ncol(kgibbs))
  colnames(schain) <- colnames(cchain) <- .multivarChainNames(snames,snames)[Kindex]
  
  snames <- otherpar$snames
  s1 <- diag(S)*0
  s2 <- r1 <- r2 <- s1
  
  message('expanding covariance chains')
  
  pbar <- txtProgressBar(min=1,max=ns,style=1)
  
  sinvPlus <-  sinvMinus <- matrix(0,S,S)   # different from zero
  
  k <- 1
  
  for(j in simIndex){
    if(REDUCT){
      Z  <- matrix(sgibbs[j,],N,r)
      ss <- .expandSigma(sigErrGibbs[j], S, Z = Z, kgibbs[j,], REDUCT = REDUCT)
      si <- .invWoodburryArma(sigErrGibbs[j], Z[kgibbs[j,],])
      cc <- .cov2Cor(ss)
      dc <- diag(sqrt(diag(ss)))
      ci <- dc%*%si%*%dc
    } else {
      ss <- .expandSigma(sgibbs[j,], S = S, REDUCT = REDUCT)
      si <- ci <- diag(1,S)
      si[notOther,notOther] <- solve(ss[notOther,notOther])
      cc <- .cov2Cor(ss)
      ci[notOther,notOther] <- solve(cc[notOther,notOther])
    }
    
    s1 <- s1 + ss
    s2 <- s2 + ss^2
    r1 <- r1 + cc
    r2 <- r2 + cc^2
    
    if(!CHAINSONLY){
      schain[k,]    <- ss[Kindex]
      cchain[k,]    <- cc[Kindex]
      if(REDUCT)kchain[k,] <- kgibbs[j,]
    }
    
    sinvPlus[si > 0]  <- sinvPlus[si > 0] + 1
    sinvMinus[si < 0] <- sinvMinus[si < 0] + 1
    
    setTxtProgressBar(pbar,k)
    k <- k + 1
  }
  diag(sinvPlus) <- diag(sinvMinus) <- 0
  sigInvPos <- which(sinvPlus > .95*length(simIndex),arr.ind=T)
  sigInvNeg <- which(sinvMinus > .95*length(simIndex),arr.ind=T)
  
  ssi <- sort( unique(c( sigInvPos[,1], sigInvNeg[,1]) ) )

  sMu  <- s1/ns
  vv   <- s2/ns - sMu^2
  vv[vv < tiny] <- tiny
  sSe  <- sqrt( vv )
  rMu  <- r1/ns
  vv   <- r2/ns - rMu^2
  vv[vv < tiny] <- tiny
  rSe  <- sqrt( vv )
  
  rownames(sMu)    <- colnames(sMu) <- snames
  rownames(sSe)    <- colnames(rSe) <- snames
  colnames(cchain) <- colnames(schain)
  
  chainList <- list(cchain = cchain, schain = schain, kchain = kchain)
  
  list(chainList = chainList, rMu = rMu, rSe = rSe, 
       sMu = sMu, sSe = sSe)
}

.expandSigma <- function(sigma, S, Z = NULL, K = NULL, REDUCT = F){
  
  if(REDUCT) return( sigma*diag(S) + tcrossprod(Z[K,]) )
  
  ss <- diag(S)
  ss[lower.tri(ss,diag=T)] <- sigma
  ss[upper.tri(ss)] <- t(ss)[upper.tri(ss)]
  ss
}

.ordTraitsFromWts <- function(yWt,ordTraits){
  
  # yWt - n by S species weights
  # ordTraits - S by p ordinal traits
  # returns n by p modal ordinal values
  
  if(!is.matrix(ordTraits))ordTraits <- matrix(ordTraits)
  
  n <- nrow(yWt)
  s <- ncol(yWt)
  
  ii <- rep(c(1:n),s)
  omat <- matrix(NA,n,ncol(ordTraits))
  
  for(j in 1:ncol(ordTraits)){
    
    PLUS <- F
    
    oj  <- ordTraits[,j]
    if(min(oj) < 0)stop('ordinal scores cannot be < 0')
    if(min(oj) == 0){
      PLUS <- T
      oj   <- oj + 1
    }
    
    rj  <- range(oj, na.rm=T)
    mm  <- matrix(0, n, rj[2] )
    jj  <- as.vector( matrix(oj, n, s, byrow=T) )
    tmp <- .byRcpp(as.vector(yWt),ii,jj,mm,mm,fun='sum')
    w0  <- which( apply(tmp,1,sum) == 0)
    
    m1  <- apply(tmp,1,which.max)
    m1  <- (rj[1]:rj[2])[m1]
    if(PLUS)m1 <- m1 - 1
    omat[,j] <- m1
    if(length(w0) > 0)omat[w0,j] <- 0
  }
  colnames(omat) <- colnames(ordTraits)
  omat
}

.incidence2Grid <- function(specs, lonLat, nx = NULL, ny = NULL, dx = NULL, 
                            dy = NULL, predGrid = NULL, effortOnly=TRUE){
  
  # must have either ngrid X 2 prediction grid, or 
  #   numbers of points nx, ny, or
  #   densities of points dx, dy
  
  ngrid <- length(predGrid)
  mapx  <- range(lonLat[,1])
  mapy  <- range(lonLat[,2])
  
  specs  <- as.character(specs)
  
  ynames <- sort(unique(specs))
  nspec  <- length(ynames)
  jj     <- match(specs,ynames)
  
  if(ngrid == 0){
    if(!is.null(dx)){
      xseq <- seq(mapx[1], mapx[2], by = dx)
      yseq <- seq(mapy[1], mapy[2], by = dy)
    } else {
      xseq <- seq(mapx[1], mapx[2], length = nx)
      yseq <- seq(mapy[1], mapy[2], length = ny)
    }
    predGrid <- as.matrix( expand.grid(lon = xseq, lat = yseq) )
    ngrid    <- nrow(predGrid)
  }
  
  ii <- RANN::nn2(predGrid, lonLat, k = 1  )$nn.idx
  mm <- matrix(0, ngrid, nspec )
  
  gridBySpec <- .byRcpp(ii*0 + 1, ii, jj, mm, mm, fun='sum')
  colnames(gridBySpec) <- ynames
  effort <- rowSums(gridBySpec)
  
  if(effortOnly){
    wk <- which(effort > 0)
    effort     <- effort[wk]
    gridBySpec <- gridBySpec[wk,]
    predGrid   <- predGrid[wk,]
  }
  list(gridBySpec = gridBySpec, predGrid = predGrid)
}


.spec2Trait <- function(pbys, sbyt, tTypes){
  
  # plotBySpec  - n by S numeric matrix
  # specByTrait - S by M data.frame
  # traitTypes  - data types for traits
  # FC can be factors that will be categorical
  
  n <- nrow(pbys)
  S <- ncol(pbys)
  M <- ncol(sbyt)
  
  ttt <- numeric(0)
  
  y2t  <- match(colnames(pbys),rownames(sbyt))  
  y2tf <- which(is.finite(y2t))
  t2y  <- match(rownames(sbyt),colnames(pbys))
  t2yf <- which(is.finite(t2y))
  
  if(is.data.frame(pbys))pbys <- as.matrix(pbys)
  
  ywt <- sweep(pbys,1,rowSums(pbys,na.rm=T),'/')
  ywt[is.na(ywt)] <- 0
  
  newTypes <- character(0)
  tmat     <- ttt <- numeric(0)
  
  ###################### neither ordinal nor factors (FC)
  
  wf   <- which(!tTypes %in% c('OC','CAT')) 

  if(length(wf) > 0){
    newTypes <- tTypes[wf]
    ttt <- sbyt[y2t,wf, drop=F]
    tmat <- ywt%*%as.matrix(sbyt[y2t,wf, drop=F])
  }
  
  ###################### ordinal classes
  
  ordNames <- which(tTypes == 'OC')
  
  if(length(ordNames) > 0){
    ordTraits <- as.matrix( round(sbyt[y2t[y2tf],ordNames],0) )
    ordCols   <- .ordTraitsFromWts(ywt,ordTraits)
    if(is.null(colnames(ordCols)))colnames(ordCols) <- colnames(ordTraits) <- 
      colnames(sbyt)[ordNames]
    ttt <- cbind(ttt, ordTraits )
    tmat <- cbind(tmat,ordCols)
    newTypes <- c(newTypes,tTypes[ordNames])
  }
  
  ##################### CAT to FC
  
  censor <- NULL
  mcol   <- ncol(tmat)
  if(is.null(mcol))mcol <- 0
  xx     <- numeric(0)
  FCgroups <- rep(0,mcol)   
  
  wf <- numeric(0)
  for(j in 1:ncol(sbyt))if(is.factor(sbyt[,j]))wf <- c(wf,j)
  
  wf <- union(wf,which(tTypes %in% 'CAT'))
  
  if(length(wf) > 0){
    
    xx <- sbyt[,wf,drop=F]
    xc <- numeric(0)
    kg <- 0
    
    for(kk in 1:length(wf)){
      
      xkk  <- xx[[kk]]            #rare type is reference
      xtab <- table(xkk)
      if(length(xtab) == 1){
        stop( paste('CAT trait _',names(xx)[kk],
                    '_ has only 1 level, need at least 2',sep='') )
      }
        
      xtab <- xtab[order(xtab)]
      xkk  <- relevel(xkk,ref=names(xtab)[1])
      cont <- contrasts(xkk,contrasts = F)
      xk   <- cont[xkk,,drop=F]
      tmp  <- ywt[,t2y]%*%xk[t2y,]
      
      if(ncol(tmp) == 2){
        mc    <- mcol + 1
        ktype <- 'CA'
        tmp <- tmp[,1,drop=F]
        gk  <- 0
        tc <- gjamCensorY( values = c(0,1), 
                           intervals = cbind( c(-Inf,0), c(1,Inf) ),
                           y = tmp)
        ttt <- cbind(ttt,xk[,1,drop=F])
        
        if(is.null(censor)){
          censor <- append(censor, list('CA' = tc$censor))
          censor$CA$columns <- mc
        } else {
          censor$CA$columns <- c(censor$CA$columns,mc)
        }
        
      } else {
        
        mc    <- ncol(tmp)
        cname <- colnames(tmp)
        cname[1] <- 'other'
        cname <- paste(colnames(xx)[kk],cname,sep='')
        colnames(tmp) <- colnames(xk) <- cname
        
        ttt   <- cbind(ttt,xk)
        ktype <- rep('FC',ncol(tmp))
        
        kg    <- kg + 1
        gk    <- rep(kg,mc)
      }
      
      mcol <- mcol + ncol(tmp)
      
      FCgroups <- c(FCgroups,gk)
      xc   <- cbind(xc,tmp)
      newTypes <- c(newTypes,ktype)

    }
    tmat <- cbind(tmat,xc)
  }
  
  colnames(tmat) <- colnames(ttt)
  
  attr(newTypes,'FCgroups') <- FCgroups

  list(plotByCWM = tmat, traitTypes = newTypes, censor = censor,
       specByTrait = ttt)
}
                

.boxplotQuant <- function( xx, ..., boxfill=NULL ){
  
  tmp <- boxplot( xx, ..., plot=F)
  ss  <- apply( xx, 2, quantile, c(.025, .158655, .5, .841345, .975) ) 
  tmp$stats <- ss
  
 pars <- list(...)
 if( 'col' %in% names(pars) )boxfill <- pars$col
  
  bxp( tmp, ..., boxfill = boxfill )
  
  tmp
}

.gjamOrd <- function( output, specLabs, col, cex, PLOT, method ){
  
  # method can be 'PCA' or 'NMDS'
  
  ematrix   <- output$parameterTables$ematrix
  ematAlpha <- output$modelSummary$ematAlpha
  
  whConZero <- output$modelSummary$whConZero
  whichZero <- output$modelSummary$whichZero
  
  y <- output$y
  S <- SO <- ncol(y)
  snames  <- colnames(y)
  
  if(is.null(col))col <- rep('black',S)
  
  other <- grep('other',colnames(y))
  
  notOther <- c(1:S)
  if(length(other) > 0){                     
    notOther   <- notOther[!notOther %in% other]
    SO         <- length(notOther)
  }
  
  plab <- c('Axis I', 'Axis II', 'Axis III')
  
  if (method == 'NMDS') {
    tmp    <- isoMDS(.cov2Dist(ematrix[notOther,notOther]), k = 3)
    eVecs  <- tmp$points
    colnames(eVecs) <- paste('NMDS',c(1:3),sep = '_')
    eValues <- lambda <- cl <- NULL
  } else {
    tmp <- eigen(ematrix[notOther,notOther])    # PCA
    eVecs   <- tmp$vectors
    eValues <- tmp$values
    lambda  <- eValues/sum(eValues)
    cl      <- cumsum(lambda)
    clab    <- paste(' (',round(100*lambda,0),'%)',sep='')
    plab    <- paste(plab, clab, sep='')
  }
  rownames(eVecs) <- snames[notOther]
  
  if(!PLOT) return( list(eVecs = eVecs, eValues = eValues) )
  
  cbord <- .getColor(col[notOther],.6)
  
  par(mfcol=c(2,2), bty='n', cex = cex, mar=c(4,4,1,1))
  
  plot(eVecs[,1],eVecs[,2],cex=1,col=cbord, bg = cbord, pch=16,
       xlab=plab[1], ylab = plab[2]) 
  abline(h=0,col=.getColor('black',.3),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.3),lwd=2,lty=2)
  
  if(length(specLabs) > 0){
    mmm <- match(specLabs,rownames(eVecs))
    text(eVecs[mmm,2],eVecs[mmm,3],specLabs,col=cbord[notOther][mmm])
  }
  
  plot(eVecs[,1],eVecs[,3],cex=1,col=cbord, bg = cbord, pch=16,
       xlab=plab[1], ylab = plab[3]) 
  abline(h=0,col=.getColor('black',.3),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.3),lwd=2,lty=2)
  
  if(length(specLabs) > 0){
    mmm <- match(specLabs,rownames(eVecs))
    text(eVecs[mmm,2],eVecs[mmm,3],specLabs,col=cbord[notOther][mmm])
  }
  
  plot(eVecs[,2],eVecs[,3],cex=1,col=cbord, bg = cbord, pch=16,
       xlab=plab[2], ylab = plab[3])
  abline(h=0,col=.getColor('black',.3),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.3),lwd=2,lty=2)
  
  if(length(specLabs) > 0){
    mmm <- match(specLabs,rownames(eVecs))
    text(eVecs[mmm,2],eVecs[mmm,3],specLabs,col=cbord[notOther][mmm])
  }
  
  if(method == 'PCA'){
    plot(cl,type='s',xlab='Rank',ylab='Proportion of variance',xlim=c(.9,S),
         ylim=c(0,1),log='x',lwd=2)
    lines(c(.9,1),c(0,cl[1]),lwd=2,type='s')
    for(j in 1:length(lambda))lines(c(j,j),c(0,cl[j]),col='grey')
    lines(cl,lwd=2,type='s')
    abline(h=1,lwd=2,col=.getColor('grey',.5),lty=2)
  }
  
  list(eVecs = eVecs, eValues = eValues)
}
                
                

