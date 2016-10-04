

gjamIIE <- function(output, xvector, MEAN = T, keepNames = NULL,
                    omitY = NULL, sdScaleX = T, sdScaleY = F){
  
  xMu    <- colMeans(output$x)
  xSd    <- apply(output$x,2,sd)
  standX <- cbind(xMu,xSd)
  colnames(standX) <- c('xmean','xsd')
  
  
  xii <- which(!names(xvector) %in% colnames(output$x))
  if(length(xii) > 0)xvector <- xvector[-xii]
  
  xii <- which(!colnames(output$x) %in% names(xvector))
  if(length(xii) > 0){
    stop('xvector is missing variables in model')
  }
  
  IIE <- .directIndirectCoeffs( snames = colnames(output$y), xvector, 
                                chains = output$chains, MEAN,
                                factorList = output$modelSummary$factorList,
                                keepNames, omitY, sdScaleY, sdScaleX,
                                standX = standX, 
                                otherpar = output$otherpar,
                                REDUCT = output$REDUCT, ng = output$modelList$ng, 
                                burnin = output$modelList$burnin)
  fit <- append(output, list('IIE' = IIE))
  fit
}