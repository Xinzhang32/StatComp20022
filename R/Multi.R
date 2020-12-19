#' @title Multi-attribute Decision Based on Information Entropy using R
#' @description Multi-attribute Decision Based on Information Entropy using R
#' @param x the data matrix which the columns represent the plans and the rows represent the variables considered
#' @return a vector of attribute values of plans which were sorted from small to large and choose the best plan
#' @examples
#' \dontrun{
#' x<-matrix(rnorm(16),4)
#' select(x)
#' }
#' @export
select<-function(x)
{
  nrow<-dim(x)[1];
  ncol<-dim(x)[2];
  R<-r<-a<-matrix(0,nrow,ncol)
  E<-numeric(ncol)
  for (i in 1:nrow)
  {
    for (j in 1:ncol)
    {
      R[i,j]<-(x[i,j]-min(x[,j]))/(max(x[,j])-min(x[,j]))
    }
  }
  
  for (i in 1:nrow)
  {
    for (j in 1:ncol)
    {
      r[i,j]<-R[i,j]/sum(R[,j])
    }
  }
  
  
  for (j in 1:ncol)
  {
    for (i in 1:nrow)
    {
      if (r[i,j]==0)
      {
        a[i,j]<-0;
      }
      else
      {
        a[i,j]<-r[i,j]*log(r[i,j])
      }
    }
    
    E[j]<--1/log(nrow)*sum(a[,j])
  }
  w<-sapply(1:ncol, function(j){
    (1-E[j])/sum(1-E[1:ncol])}
  )
  z<-sapply(1:nrow,function(i) {
    sum(r[i,]*w)
  })
  return(c(sort(z),max(sort(z)),which(z==max(sort(z)))))
}
