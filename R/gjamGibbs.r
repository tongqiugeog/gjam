
gjamGibbs <- function(formula, xdata, ydata, modelList){
  
  # xdata        - n by Q design data.frame
  # ydata        - n by S response, continuous, discrete, categorical
  
  .gibbsLoop( formula, xdata, ydata, modelList)
}
  