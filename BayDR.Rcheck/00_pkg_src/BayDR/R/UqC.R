#' Decompose the confounder matrix
#'
#' This function decomposes the confounder matrix.
#' 
#' 
#' @param Y binary repsonse vector length of n
#' @param X binary exposure vector same length with Y
#' @param C binary confounder matrix n*r, r the number of counfouders
#' @param Dat can be missing if given Y,X,C
#' @return a list UqC of the unique levels of C, N..Number, N..Number0, N..Number1 
#' @export
#' @importFrom plyr dlply
#' @examples
#' \dontrun{
#' data(sampleDat)
#' Y=sample.dataset$Y
#' X=sample.dataset$X
#' C=sample.dataset$C
#' dcom(Y,X,C)
#' }

dcom<- function(C,X,Y,Dat){
  
  if(missing(Dat)) Dat=cbind(C,X,Y)
  
  Dat <- as.data.frame(Dat)
  no.confounder = ncol(Dat)-2
  names(Dat) <- c(paste("C",1:no.confounder,sep=""),"X","Y")
  Dat$X <- factor(Dat$X)
  Dat$Y <- factor(Dat$Y)
  temp.list <- dlply(Dat,paste("C",1:no.confounder,sep=""),
                     function(x) {
                       temp.table <- table(x[,(no.confounder+1):(no.confounder+2)])
                       as.numeric(temp.table)
                     }
  )
  UqC <- attr(temp.list,"split_labels")
  N..Number <- as.data.frame(do.call(rbind,lapply(temp.list,function(x) c(x[4],x[4]+x[2],x[3],x[3]+x[1],sum(x)))))
  
  colnames(N..Number) <- c('C11','C1.','C01','C0.','Number') # matrix for n_cxx 
  rownames(N..Number) <- 1:nrow(N..Number)
  N..Number0 <- N..Number[,3:4];N..Number1 <- N..Number[,1:2] # for C=0 & for C=1
  return(
    list(UqC=UqC,N..Number=N..Number,N..Number0=N..Number0,N..Number1=N..Number1)
  )
  
}