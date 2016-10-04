

gjamPoints2Grid <- function(specs, xy, nxy = NULL, dxy = NULL, 
                            predGrid = NULL, effortOnly = TRUE){
  
  dx <- dy <- nx <- ny <- NULL
  
  if(length(nxy) == 1) nx <- ny <- nxy
  if(length(dxy) == 1) dx <- dy <- dxy
  if(is.null(nxy) & is.null(dxy) & is.null(predGrid)){
    stop('must supply nxy or dxy or predGrid')
  }
  if(length(nxy) == 2) nx <- nxy[1]; ny <- nxy[2]
  if(length(dxy) == 2) dx <- dxy[1]; dy <- dxy[2]
  
  if(length(specs) != nrow(xy))stop('specs must have length = nrow(xy)')
  
  mr <- apply( apply(xy,2,range), 2, diff )
  
  if(!is.null(dxy)){
    if(mr[1] < (3*dx) | mr[2] < (3*dy))stop('dxy too small')
  }
  .incidence2Grid(specs, lonLat = xy, nx, ny, dx, dy, predGrid, effortOnly)
}
