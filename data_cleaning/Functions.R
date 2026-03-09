##############
#FUNCTIONS
##############

#' Compute Box-Cox Transformation for a given lambda
#'
#' @param x Vector of strictly positive data.
#' @param lambda Chosen lambda based on profile likelihood plot
#'
#' @return Box-Cox transformed x.
#' 
boxcox_transform <- function(x, lambda){
  
  if(sum(x <= 0) > 0){
    stop("Box-Cox Transform only appropriate for strictly positive data.
         Some observations <= 0.")
  }
  if(lambda == 0){
    return(log(x))
  } else{
    return((x^lambda - 1) / lambda)
  }
  
}


#' Calculate standard errors
#' @param x Vector of data.

calcSE<-function(x){
  x <- x[is.na(x)==F] 
  sd(x)/sqrt(length(x))
}