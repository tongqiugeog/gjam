
gjamTrimY <- function(y, minObs = 2, maxCols = NULL, OTHER = T){  
    
    # minObs    - minimum no. of non-zero values in a column of y
    # maxCols   - number of columns to retain, those with highest values
    # if(OTHER) sum of rare are returned in 'other' column
    
    y      <- as.matrix(y)
    nc     <- ncol(y)
    mnames <- colnames(y)
    
    io <- y
    io[io > 0] <- 1
    
    csum <- colSums(io, na.rm=T)
    ww <- which(csum > minObs)
    
    if(!is.null(maxCols)){
      ww <-  ww[ order(csum[ww],decreasing=T) ]
      ww <-ww[1:maxCols]
    }
    
    out    <- y[,ww]
    mnames <- mnames[ww]
    
    if(OTHER){
      other <- rowSums(y[,-ww],na.rm=T)
      out <- cbind(out,other)
      mnames <- c(mnames,'other')
    }
    
    if(!is.matrix(out)){
      out <- matrix(out,ncol=1)
      colnames(out) <- mnames
    }
    csum <- csum[ww]
    
    list(y = out, colIndex = ww, nobs = csum)
  }
  