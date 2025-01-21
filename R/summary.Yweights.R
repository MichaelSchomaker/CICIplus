summary.Yweights<-function(object,verbose=TRUE,...){

x <- object
time.points <- length(x[[1]])
n.int <- ncol(x[[1]][[1]]) 
n.c <- length(x)
  
prop.not.one <- function(vec){mean(vec!=1)}
not.one      <- function(mat){apply(mat,2,prop.not.one)}
mymin        <- function(mat){apply(mat,2,min)}
mymax        <- function(mat){apply(mat,2,max)}

tab <- matrix(NA,ncol=n.int,nrow=3,dimnames = list(c("% weights not 1","min(weight)","max(weight)"),
                                                   colnames(x[[1]])))
results <- rep(list(rep(list(tab),time.points)),n.c); names(results)<-names(x)

for(i in 1:n.c){
proportion.not.one <- lapply(x[[i]], not.one)
savemin <- lapply(x[[i]],mymin)
savemax <- lapply(x[[i]],mymax)
names(results[[i]]) <- paste("Time point", 1:time.points)
  for(j in 1:time.points){
    results[[i]][[j]][1,] <- proportion.not.one[[j]]*100
    results[[i]][[j]][2,] <- savemin[[j]]
    results[[i]][[j]][3,] <- savemax[[j]]
    results[[i]][[j]]     <- round(results[[i]][[j]], digits=2)
  }
}
if(verbose==TRUE){print(results)}
return(invisible(results))
}
  