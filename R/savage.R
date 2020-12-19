#' @title Choose the best plan based on Savage using R
#' @description Choose the best plan based on Savage using R
#' @param x the data matrix which the columns represent the plans and the rows represent the variables considered 
#' @return the best solution based on Savage criteria and its compromise standard return value
#' @examples
#' \dontrun{
#' x<-matrix(c(30,20,10,12,15,9,-6,-2,5),3)
#' savage.select(x)
#' }
#' @export
savage.select<-function(x)
{
  nrow<-dim(x)[1]
  ncol<-dim(x)[2]
  savage<-sapply(1:ncol,function(j) {
    sapply(1:nrow,function(i) {
      max(x[,j])-x[i,j]
    })  
  })
  choose<-sapply(1:nrow,function(i){
    max(savage[i,]) })
  return(c(min(choose),which(choose==min(choose))))
}
