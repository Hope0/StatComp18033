#' @title quicksort and count using R.
#' @description quicksort is a efficiency but not a robust algorithm, the function mainly shows the exchange and comparison times of quicksort
#' @param x the vector to be sorted
#' @param is.count logical; if TRUE, return the times of comparisons and exchanges
#' @param decreasing  logical; the sort be increasing or decreasing
#' @return if is.count=FALSE, the output is the result as a vector, else if is.count=TRUE, the output is a list with result ,compare and exchange.
#' @examples
#' \dontrun{
#' X<-rnorm(10,0,1)
#' quicksort(X)
#' quicksort(X,is.count=TRUE)
#' }
#' @export

quicksort <- function(x,is.count=FALSE,decreasing = FALSE){
  quicksort1 <- function(x){
    if(length(x)<2) return(x)
    r <- x[1]
    xsn <- x<r
    xs <- x[xsn ]
    xl <- x[!xsn]
    xl <- xl[-1]
    return(c(quicksort(xs),r,quicksort(xl)))
  }
  quicksortcs <- function(X){
    compare <- X[1]
    exchange <- X[2]
    x<- X[c(-1,-2)]
    if(length(x)<2) return(c(compare,exchange,x))
    r <- x[1]
    xsn <- x<r
    xs <- x[xsn ]
    xl <- x[!xsn]
    xl <- xl[-1]
    left <- quicksortcs(c(0,0,xs))
    right <- quicksortcs(c(0,0,xl))
    compare =compare +length(x)-1 + left[1] +right[1]
    exchange =exchange +sum(xsn)+left[2] +right[2]
    return(c(compare,exchange,left[c(-1,-2)],r,right[c(-1,-2)]))
  }
  x<-as.vector(x)
  if(decreasing ==FALSE){
    if(is.count==FALSE) return(quicksort1(x))
    if(is.count==TRUE){
      ans <- quicksortcs(c(0,0,x))
      return(list(result=ans[c(-1,-2)],compare=ans[1],exchange=ans[2]))
    }
  }
  if(decreasing ==TRUE){
    x <- -x
    if(is.count==FALSE) return(-quicksort1(x))
    if(is.count==TRUE){
      ans <- quicksortcs(c(0,0,x))
      return(list(result=-ans[c(-1,-2)],compare=ans[1],exchange=ans[2]))
    }
  }
  print("error")
}
