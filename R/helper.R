linear <- function(nom,den,c){(apply(nom > c,2,as.numeric))*1 +  (apply(nom <= c,2,as.numeric))*(nom/den)}

multiResultClass <- function(result1=NULL,result2=NULL)
{
  me <- list(
    result1 = result1,
    result2 = result2
  )
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

# 1) binning
binning <- function(form.n,form.d,X,Anodes,abar){
  #
  if(length(abar)<2){stop("abar needs to have >= 2 values to apply binning")}
  contin_var <- apply(subset(X, select=Anodes), 2, function(var) length(unique(var)))
  if(any(contin_var<10)){warning("Some of your intervention variables have less than 10 unique values. Are you sure A is continuous?")}
  #
  fitted.n <- fitted.d <- rep(list(NULL),length(form.n))
  g.d <- g.n <- rep(list(matrix(NA,nrow=nrow(X),ncol=length(abar),dimnames=list(NULL,paste(abar)))),
                    length(form.n))
  cuts <- rep(NA,length(abar)-1)
  for(i in 1:(length(cuts))){cuts[i] <- (abar[i] + abar[i+1])/2}  
  #cuts <- c(-Inf,cuts,Inf)
  cuts <- c(cuts[min(order(cuts))]-mean(diff(abar)),cuts,cuts[max(order(cuts))]+mean(diff(abar)))
  #
  for(i in 1:length(form.n)){
    for(j in 1:length(abar)){
      cutX <- X
      cutX[,Anodes] <- apply((subset(cutX, select=Anodes)>cuts[j]) & (subset(cutX, select=Anodes) < cuts[j+1]),2,as.numeric)
      fitted.n[[i]] <- mgcv::gam(as.formula(form.n[i]),data=cutX,family="binomial")
      fitted.d[[i]] <- mgcv::gam(as.formula(form.d[i]),data=cutX,family="binomial")
      g.n[[i]][,j] <- dbinom(1,size=1,predict(fitted.n[[i]], type="response", newdata=X))
      g.d[[i]][,j] <- dbinom(1,size=1,predict(fitted.d[[i]], type="response", newdata=X))
    }
  }
  #
  return(list(g.n,g.d))
  #
}

# 2) dnorm, dpois, dbinom

parametric <- function(form.n,form.d,X,Anodes,abar){
  fams <- CICI:::assign.family(as.data.frame(X[,Anodes]))
  if(any(substr(fams,1,4)=="mult")){stop("multinomial intervention currently not supported")}
  fitted.n <- fitted.d <- rep(list(NULL),length(fams))
  g.d <- g.n <- rep(list(matrix(NA,nrow=nrow(X),ncol=length(abar),dimnames=list(NULL,paste(abar)))),
                         length(fams))
  for(i in 1:length(form.n)){
  fitted.n[[i]] <- mgcv::gam(as.formula(form.n[i]),data=X,family=fams[i])
  fitted.d[[i]] <- mgcv::gam(as.formula(form.d[i]),data=X,family=fams[i])
    for(j in 1:length(abar)){
      if(fams[i]=="gaussian"){
      g.n[[i]][,j] <- dnorm(abar[j],mean=predict(fitted.n[[i]], newdata=X),sd=sqrt(fitted.n[[i]]$sig2))
      g.d[[i]][,j] <- dnorm(abar[j],mean=predict(fitted.d[[i]], newdata=X),sd=sqrt(fitted.d[[i]]$sig2))
      }else{
        if(fams[i]=="poisson"){
          g.n[[i]][,j] <- dpois(abar[j],lambda=predict(fitted.n[[i]], type="response", newdata=X))
          g.d[[i]][,j] <- dpois(abar[j],lambda=predict(fitted.d[[i]], type="response", newdata=X))
        }else{
          g.n[[i]][,j] <- dbinom(abar[j],size=1,prob=predict(fitted.n[[i]], type="response", newdata=X))
          g.d[[i]][,j] <- dbinom(abar[j],size=1,prob=predict(fitted.d[[i]], type="response", newdata=X))
            }
      }
  }}
  return(list(g.n,g.d))
}

# 3) haldensify

hal_density <- function(form.n,form.d,X,Anodes,abar,...){
  #
  if(length(abar)<2){stop("abar needs to have >= 2 values to apply haldensify")}
  contin_var <- apply(subset(X, select=Anodes), 2, function(var) length(unique(var)))
  if(any(contin_var<10)){warning("Some of your intervention variables have less than 10 unique values. Are you sure A is continuous?")}
  #
  fitted.n <- fitted.d <- rep(list(NULL),length(form.n))
  g.d <- g.n <- rep(list(matrix(NA,nrow=nrow(X),ncol=length(abar),dimnames=list(NULL,paste(abar)))),
                    length(form.n))
  #
  decide_knots_1 = ifelse(ncol(X)>=20,c(25, 10),c(25, 10, 5))
  decide_knots_0 = ifelse(ncol(X)>=20,c(50, 25),c(50, 25, 10))
  for(i in 1:length(form.n)){
        #
        fitted.n[[i]] <- haldensify::haldensify(
          A = as.numeric(unlist(c(subset(X,select=gsub(" ", "",strsplit(form.n[i],"~")[[1]][1]))))),
          W = subset(X,select=gsub(" ", "",strsplit(strsplit(form.n[i],"~")[[1]][2],"[+]")[[1]])),
          n_bins = 100,
          grid_type = c("equal_range","equal_mass"),
          lambda_seq = exp(seq(-0.1, -10, length = 100)),
          # 
          num_knots = hal9001:::num_knots_generator(
            max_degree = ifelse(ncol(X)>=20,2,3),
            smoothness_order = 1,
            base_num_knots_0 = decide_knots_0,
            base_num_knots_1 = decide_knots_1
          )
        )
        #
        tc <- gsub(" ", "",strsplit(strsplit(form.d[i],"~")[[1]][2],"[+]")[[1]])
        if(tc[1]=="1" & length(tc)==1){W.d<-rep(1,nrow(X))}else{
          W.d<- subset(X,select=gsub(" ", "",strsplit(strsplit(form.d[i],"~")[[1]][2],"[+]")[[1]])[-1])}
        
        fitted.d[[i]] <-haldensify::haldensify(
          A = as.numeric(unlist(c(subset(X,select=gsub(" ", "",strsplit(form.d[i],"~")[[1]][1]))))),
          W = W.d,
          n_bins = 100,
          grid_type = c("equal_range","equal_mass"),
          lambda_seq = exp(seq(-0.1, -10, length = 100)),
          # 
          num_knots = hal9001:::num_knots_generator(
            max_degree = ifelse(ncol(X)>=20,2,3),
            smoothness_order = 1,
            base_num_knots_0 = decide_knots_0,
            base_num_knots_1 = decide_knots_1
          )
        )
        #
  for(j in 1:length(abar)){
      g.n[[i]][,j] <- haldensify:::predict.haldensify(fitted.n[[i]], new_A = rep(abar[j],nrow(X)),
                                             new_W = subset(X,select=gsub(" ", "",strsplit(strsplit(form.n[i],"~")[[1]][2],"[+]")[[1]])), trim=FALSE, ...)
      g.d[[i]][,j] <- haldensify:::predict.haldensify(fitted.d[[i]], new_A = rep(abar[j],nrow(X)),
                              new_W = subset(X,select=gsub(" ", "",strsplit(strsplit(form.n[i],"~")[[1]][2],"[+]")[[1]])), trim=FALSE, ...)
    }
  }
  #
  return(list(g.n,g.d))
  #
}

###

find.Qs <- function(dat, Y, A, L, C){
  position.AC <- which(colnames(dat)%in%c(A, C))
  position.LY <- which(colnames(dat)%in%c(L, Y))
  position.Y  <- which(colnames(dat)%in%c(Y))
  Q.list <- rep(list(NA), length(position.Y))
  maxA.list <- rep(list(NA), length(position.Y)) 
  block.index <- rep(NA,length(position.LY))
  block.index.AC <- rep(NA,length(position.AC))
  block.index[length(block.index)]<-1; block.index.AC[length(block.index.AC)]<-1
  if(length(position.LY)>1){for(j in length(position.LY):2){
    if(position.LY[j]-position.LY[j-1]>1){block.index[j-1]<-block.index[j]+1}else{block.index[j-1]<-block.index[j]}
  }}
  Qs <- colnames(dat)[sort(unlist(lapply(split(position.LY,block.index),min)))]
  if(length(position.AC)>1){for(j in length(position.AC):2){
    if(position.AC[j]-position.AC[j-1]>1){block.index.AC[j-1]<-block.index.AC[j]+1}else{block.index.AC[j-1]<-block.index.AC[j]}
  }}
  max.ACblocks <- colnames(dat)[sort(unlist(lapply(split(position.AC,block.index.AC),max)))]
  # 
  for(i in 1:length(Q.list)){
  lastY <- colnames(dat)[position.Y[i]]
  dat2 <- dat[,1:position.Y[i]]
  mod.position.AC <- position.AC[which(position.AC<position.Y[i])]
  all.relevant.Qs <- Qs[Qs%in%colnames(dat2)]
  unallowed.Qs <- colnames(dat2)[c(1:min(mod.position.AC),max(mod.position.AC):ncol(dat2))]
  if(any(all.relevant.Qs%in%unallowed.Qs)){Qsi <- all.relevant.Qs[!all.relevant.Qs%in%unallowed.Qs]}else{Qsi<-all.relevant.Qs}
  Q.list[[i]] <- c(Qsi,lastY)
  maxA.list[[i]] <- rep(NA, length(Q.list[[i]]))
    for(k in 1:length(maxA.list[[i]])){
        maxA.list[[i]][k] <- colnames(dat)[max(position.AC[which(position.AC<which(colnames(dat2)%in%Q.list[[i]][k]))])]}
  }
  #
  AC.block <- AC.block.names <-  split(position.AC,block.index.AC)
  for(i in 1:length(AC.block.names)){AC.block.names[[i]] <- colnames(dat)[AC.block[[i]]] }
  #
  return(list(Qs.for.each.Y=Q.list, maxAC.for.each.Q=maxA.list, AC.block.names=AC.block.names))
} 

Multiply <- function(lists){
  if(length(lists)==1){output<-lists}
  if(length(lists)>=1){
    output <- Reduce("*",lists)
  }
  return(output)
}


