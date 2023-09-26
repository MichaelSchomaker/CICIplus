calc.weights <- function(X,Anodes=NULL,Ynodes=NULL,Lnodes=NULL,Cnodes=NULL,
                         abar=NULL,times=length(Anodes),c=0.01,
                         screen = FALSE, 
                         survival = FALSE, eps = 1e-10, zero=0,
                         d.method=c("binning","parametric","hal_density"),
                         z.method=c("density","eps"),
                         w.function="linear",
                         for.sgf=TRUE,
                         verbose=TRUE, ...
                         ){
# prepare  
d.method <- match.arg(d.method)
z.method <- match.arg(z.method)
if(is.null(abar)){stop("please provide values for abar")}
if(min(which(colnames(X)%in%Ynodes))<min(which(colnames(X)%in%Anodes))){stop("Ynodes can not occur before Anodes.\n Likely you specified pre-intervention variables as Ynode?")}
if(any(c(Ynodes,Lnodes,Anodes,Cnodes)%in%colnames(X)==FALSE)){stop(paste("You have specified the following variable name(s) which are not part of the data:",
                                                                         paste(c(Ynodes,Lnodes,Anodes,Cnodes)[c(Ynodes,Lnodes,Anodes,Cnodes)%in%colnames(X)==FALSE],collapse=" ") ,"\n"))}
if(max(which(colnames(X)%in%Anodes))==ncol(X)){stop("Anodes/Cnodes should not be in last column")}

# construct model formulas
Aform.n <-  CICI::make.model.formulas(X=X, Anodes=Anodes, Cnodes=Cnodes, Ynodes= Ynodes, survival=survival) 
if(screen==T){Aform.n <- CICI::model.formulas.update(Aform.n$model.names, X, pw=verbose) # ... problem with gdf() below
              Aform.n <- Aform.d <- Aform.n$Anames 
}else{Aform.n <- Aform.d <- Aform.n$model.names$Anames}

Aform.d2 <- rep(NULL,length(Aform.d))
  for(j in 1:length(Aform.d)){
    if(j>1){cond1 <- (which(colnames(X)%in%c(Lnodes,Ynodes)) >
    which(colnames(X)%in%gsub(" ", "",unlist(strsplit(Aform.d[j-1],"~"))[1])))}else{cond1<-rep(TRUE,length(c(Lnodes,Ynodes)))} 
    cond2 <- (which(colnames(X)%in%c(Lnodes,Ynodes)) <
    which(colnames(X)%in%gsub(" ", "",unlist(strsplit(Aform.d[j],"~"))[1])))
    rem <- colnames(X)[colnames(X)%in%c(Lnodes,Ynodes)][which(cond1 & cond2)]
    if(length(rem)>0){nf <- gsub(" ", "",Aform.d[j]); nf <- gsub("~","~1+",nf)
      for(k in 1:length(rem)){
        nf <- gsub(paste0("+",rem[k]),"",nf,fixed=T) 
      }}else{if(j==1){nf<-gsub(" ", "",Aform.d[j]); nf <- gsub("~","~1+",nf)}else{
                      nf<-paste0(unlist(strsplit(Aform.d[j],"~"))[1],"~",unlist(strsplit(Aform.d2[j-1],"~"))[2])}
             }
    Aform.d2[j] <- nf
    }

# estimate conditional densities
gdf <- get(d.method)
g.preds <- gdf(Aform.n,Aform.d2,X=X,Anodes=Anodes,abar=abar,...)
g.preds.n <- g.preds[[1]]; g.preds.d <- g.preds[[2]]

# replace denominator if zero
if(any(unlist(lapply(g.preds.d,function(x){any(x<=zero, na.rm=T)})))){
  if(verbose==TRUE){message("Note: zeros in denominator. Will be dealt with as specified in `z.method'.")}
  if(z.method=="density"){
    problem <- TRUE; k=1; last.resort <- NULL
    while(problem==TRUE){
    which.times <- unlist(lapply(g.preds.d,function(x){any(x<=zero)}))
    which.columns <- lapply(lapply(g.preds.d,apply,2,function(x){x<=zero}),apply,2,any)
    for(j in 1:length(g.preds.d)){
      if(which.times[j]==TRUE){
        if(j>k){Aform.d2[j] <- paste0(strsplit(Aform.d2[j],"~")[[1]][1],"~",strsplit(Aform.d2[j-k],"~")[[1]][2])}
        if(j<=k){last.resort <- c(last.resort,list(c(j,k)))}
      }
    }
    g.preds.d.update <- gdf(Aform.n,Aform.d2,X=X,Anodes=Anodes,abar=abar,...)[[2]]
    if(is.null(last.resort)==FALSE){
      for(q in 1:length(last.resort)){
      g.preds.d.update[[last.resort[[q]][1] ]][,last.resort[[q]][2]]  <- rep(1,nrow(g.preds.d[[1]]))
      }}
    for(l in 1:length(g.preds.d)){for(m in 1:ncol(g.preds.d[[1]])){if(which.columns[[l]][m]==TRUE){
      g.preds.d[[l]][,m] <- g.preds.d.update[[l]][,m]
      if(verbose==TRUE){message(paste0("Denominator zero at time ",l,", int=", colnames(g.preds.d[[1]])[m], ", try #",k))}
    }}}
    check.if.problem <- unlist(lapply(g.preds.d,function(x){any(x<=zero)}))
    if(any(check.if.problem)==TRUE){k<-k+1}else{problem <- FALSE}
    }
  }else{
    if(z.method=="eps"){g.preds.d <- lapply(g.preds.d,function(x){x[x<=zero]<-eps; return(x)})}
    }
  }

# calculate weights
wf <- get(w.function)
w.list <- rep(list(NA),length(c)); names(w.list)<-paste(c)
  for(k in 1:length(c)){
    for(j in 0:(times-1)){
    assign(paste0("w",j), wf(nom = g.preds.n[[j+1]], den = g.preds.d[[j+1]], c=c[k]))
    }
  w.list[[k]] <- mget(paste0("w",0:(times-1)))
  }

if(for.sgf==TRUE){
  loop.Q    <- find.Qs(dat=X, L=Lnodes, Y=Ynodes, A=Anodes, C=Cnodes) 
  max.Qs    <- max(unlist(lapply(loop.Q[[1]],length)))
  if(length(w.list[[1]])!=max.Qs){
    if(verbose==TRUE){message(paste0("Weights within blocks of A/C nodes are multiplied, where necessary.\n"))}
    if(length(loop.Q$AC.block.names)!=max.Qs){stop("Contact maintainer; something does not make sense.")}
    if(length(unlist(loop.Q$AC.block.names))!=length(w.list[[1]])){stop("Contact maintainer; something does not make sense.")}
    w.binding <- rep(NA,length(w.list[[1]]))
    ll <- lapply(loop.Q$AC.block.names, length)
    k <- length(w.binding); block <- 1; group<- length(loop.Q$AC.block.names)
    while(k>=1){
      w.binding[k:(k+1-ll[[block]])]<-group
      k<-k-ll[[block]] ; block <- block+1; group<-group-1
    }
    for(m in 1:length(w.list)){w.list[[m]] <- lapply(split(w.list[[m]],w.binding),Multiply) }
  }
}

weight<-w.list


# return
class(weight) <- "Yweights"
return(weight)
}