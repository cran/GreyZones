#' @export
detectGreyZones.default <- function(table ){

    R <- nrow(table)
    C <- ncol(table)
    if (R < 2){
      stop("The agreement table must have atleast two rows!")
    }
    if (C < 2){
      stop("The agreement table must have atleast two columns!")
    }
    if (R != C){
      stop("The agreement table must be a square contingency table!")
    }
    if (any(table < 0)){
      stop("The counts in the agreement table must be nonnegative!")
    }
    
    res <- detectGreyZones.main(table = table)
   
    res$call <- match.call()
    class(res) <- c("detectGreyZones")
    res
    
  }