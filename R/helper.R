gal_ga  <- function(nom,den,c){(apply(nom > c,2,as.numeric))*1 +  (apply(nom <= c,2,as.numeric))*(nom/den)}
gal_ga2 <- function(nom,den,c){(apply((nom/den) > c,2,as.numeric))*1 +  (apply((nom/den) <= c,2,as.numeric))*(nom/den)}

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
      pX <- cutX; pX[,Anodes] <- 1
      g.n[[i]][,j] <- dbinom(1,size=1,predict(fitted.n[[i]], type="response", newdata=pX))
      g.d[[i]][,j] <- dbinom(1,size=1,predict(fitted.d[[i]], type="response", newdata=pX))
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
      XA <- X; XA[,Anodes] <- abar[j] 
      if(fams[i]=="gaussian"){
        g.n[[i]][,j] <- dnorm(abar[j],mean=predict(fitted.n[[i]], newdata=XA),sd=sqrt(fitted.n[[i]]$sig2))
        g.d[[i]][,j] <- dnorm(abar[j],mean=predict(fitted.d[[i]], newdata=XA),sd=sqrt(fitted.d[[i]]$sig2))
      }else{
        if(fams[i]=="poisson"){
          g.n[[i]][,j] <- dpois(abar[j],lambda=predict(fitted.n[[i]], type="response", newdata=XA))
          g.d[[i]][,j] <- dpois(abar[j],lambda=predict(fitted.d[[i]], type="response", newdata=XA))
        }else{
          g.n[[i]][,j] <- dbinom(abar[j],size=1,prob=predict(fitted.n[[i]], type="response", newdata=XA))
          g.d[[i]][,j] <- dbinom(abar[j],size=1,prob=predict(fitted.d[[i]], type="response", newdata=XA))
            }
      }
  }}
  return(list(g.n,g.d))
}

# 3) haldensify

hal_density <- function(form.n,form.d,X,Anodes,abar,
                        n_bins=max(10, sqrt(nrow(X))),
                        lambda_seq = exp(seq(-0.1, -10, length = 100)),
                        grid_type=c("equal_range","equal_mass"),
                        max_degree = ifelse(ncol(X)>=20,2,3),
                        smoothness_orders = 1,
                        hal.verbose=FALSE,
                        runtime=c("fast","very_fast","fairly_fast","somewhat_fast","regular","not_specified"),
                        ...){
  #
  if(length(abar)<2){stop("abar needs to have >= 2 values to apply haldensify")}
  contin_var <- apply(subset(X, select=Anodes), 2, function(var) length(unique(var)))
  if(any(contin_var<10)){warning("Some of your intervention variables have less than 10 unique values. Are you sure A is continuous?")}
  runtime <- match.arg(runtime)
  very_fast        <- list(c(50,25),c(50,25,10),c(25,10),c(25, 10, 5))
  fast             <- list(c(100,50),c(100,50,25),c(40,15),c(40, 15, 10))
  fairly_fast      <- list(c(200,100),c(200,100,50),c(50,25),c(50, 25, 15)) 
  somewhat_fast    <- list(c(400,200),c(400,200,100),c(100,75),c(100, 75, 50))
  regular          <- list(c(500,200),c(500,200,50),c(200,100),c(200, 100, 50))
  sel_rt <- get(runtime)
  #
  fitted.n <- fitted.d <- rep(list(NULL),length(form.n))
  g.d <- g.n <- rep(list(matrix(NA,nrow=nrow(X),ncol=length(abar),dimnames=list(NULL,paste(abar)))),
                    length(form.n))
  #
  for(i in 1:length(form.n)){
        # numerator
        W.n <- subset(X,select=gsub(" ", "",strsplit(strsplit(form.n[i],"~")[[1]][2],"[+]")[[1]]))
        nknots <- num_knots_generator(
          max_degree = max_degree,
          smoothness_orders  = smoothness_orders,
          base_num_knots_0 = ifelse(ncol(W.n)>=20,sel_rt[[1]],sel_rt[[2]]),
          base_num_knots_1 = ifelse(ncol(W.n)>=20,sel_rt[[3]],sel_rt[[4]])
        )
        fitted.n[[i]] <- haldensify::haldensify(
          A = as.numeric(unlist(c(subset(X,select=gsub(" ", "",strsplit(form.n[i],"~")[[1]][1]))))),
          W = W.n,
          n_bins = n_bins,
          grid_type = grid_type,
          lambda_seq = lambda_seq,
          num_knots = nknots, ...
        )
        if(hal.verbose==TRUE){message(paste("Number of knots for numerator density:", paste(nknots, collapse=","), "(", form.n[i],")"))}
        # denominator
        tc <- gsub(" ", "",strsplit(strsplit(form.d[i],"~")[[1]][2],"[+]")[[1]])
        if(tc[1]=="1" & length(tc)==1){W.d<-matrix(1,nrow(X),ncol=1)}else{
        W.d<- subset(X,select=gsub(" ", "",strsplit(strsplit(form.d[i],"~")[[1]][2],"[+]")[[1]])[-1])}
        nknots2 <- num_knots_generator(
          max_degree = max_degree,
          smoothness_orders  = smoothness_orders,
          base_num_knots_0 = ifelse(ncol(W.d)>=20,sel_rt[[1]],sel_rt[[2]]),
          base_num_knots_1 = ifelse(ncol(W.d)>=20,sel_rt[[3]],sel_rt[[4]])
        )
        fitted.d[[i]] <-haldensify::haldensify(
          A = as.numeric(unlist(c(subset(X,select=gsub(" ", "",strsplit(form.d[i],"~")[[1]][1]))))),
          W = W.d,
          n_bins = n_bins,
          grid_type = grid_type,
          lambda_seq = lambda_seq,
          num_knots = nknots2,...
        )
        if(hal.verbose==TRUE){message(paste("Number of knots for denominator density:", paste(nknots2, collapse=","), "(", form.d[i],")"))}
        #
  for(j in 1:length(abar)){
      WA.n <- subset(X,select=gsub(" ", "",strsplit(strsplit(form.n[i],"~")[[1]][2],"[+]")[[1]]))
      if(length(strsplit(form.d[i],"~")[[1]][2])>1){
      WA.d <- subset(X,select=gsub(" ", "",strsplit(strsplit(form.d[i],"~")[[1]][2],"[+]")[[1]][-1]))}else{
      WA.d <- matrix(1,nrow(X),ncol=1)  
      }
      sel.n <- colnames(WA.n)%in%Anodes; WA.n[,sel.n] <- abar[j]
      sel.d <- colnames(WA.d)%in%Anodes; WA.d[,sel.d] <- abar[j]
      g.n[[i]][,j] <- haldensify:::predict.haldensify(fitted.n[[i]], new_A = rep(abar[j],nrow(X)), new_W = WA.n, trim=FALSE, ...)
      g.d[[i]][,j] <- haldensify:::predict.haldensify(fitted.d[[i]], new_A = rep(abar[j],nrow(X)), new_W = WA.d, trim=FALSE, ...)
  }
  }
  #
  return(list(g.n,g.d))
  #
}

# 4) SL density
hazardbinning <- function(form.n,form.d,X,Anodes,abar,SL.library=NULL){
  #
  if(length(abar)<2){stop("abar needs to have >= 2 values to apply binning")}
  contin_var <- apply(subset(X, select=Anodes), 2, function(var) length(unique(var)))
  if(any(contin_var<10)){warning("Some of your intervention variables have less than 10 unique values. Are you sure A is continuous?")}
  #
  g.d <- g.n <- rep(list(matrix(NA,nrow=nrow(X),ncol=length(abar),dimnames=list(NULL,paste(abar)))),
                    length(form.n))
  #
  cuts <- (head(abar, -1) + tail(abar, -1)) / 2
  cuts <- c(cuts[1] - mean(diff(abar)), cuts, cuts[length(cuts)] + mean(diff(abar)))
  #
  result_list <- lapply(seq_along(form.n), function(i) {
    A <- X[, Anodes[i]]
    fml.n <- as.formula(sub("~\\s*", "~ s(bin_id) + ", form.n[[i]]));vars.n <- all.vars(fml.n[[3]]);y.n <- all.vars(fml.n[[2]])
    fml.d <- as.formula(sub("~\\s*", "~ s(bin_id) + ", form.d[[i]]));vars.d <- all.vars(fml.d[[3]]);y.d <- all.vars(fml.d[[2]])
    W.n <- X[, vars.n[-1], drop = FALSE]
    W.d <- X[, vars.d[-1], drop = FALSE]
    y_bins0 <- findInterval(A, cuts, rightmost.closed = TRUE)
    #
    dat.n <- cbind(data.frame(y = y_bins0, status = 1), W.n)
    dat.d <- cbind(data.frame(y = y_bins0, status = 1), W.d)
    #
    long_data_s.n <- as_ped(survival::Surv(y, status) ~ ., data = dat.n, cut = c(0:length(cuts)))
    long_data_s.d <- as_ped(survival::Surv(y, status) ~ ., data = dat.d, cut = c(0:length(cuts)))
    names(long_data_s.n)[names(long_data_s.n) == "ped_status"] <- Anodes[i]
    names(long_data_s.d)[names(long_data_s.d) == "ped_status"] <- Anodes[i]
    
    
    if (!is.null(SL.library)) {
      fit_n <- SuperLearner(
        Y = long_data_s.n[, y.n],
        X = long_data_s.n[, vars.n, drop = FALSE],
        id = long_data_s.n$id,
        family = binomial(),
        verbose = FALSE,
        SL.library = SL.library
      )
      fit_d <- SuperLearner(
        Y = long_data_s.d[, y.d],
        X = long_data_s.d[, vars.d, drop = FALSE],
        id = long_data_s.d$id,
        family = binomial(),
        verbose = FALSE,
        SL.library = SL.library
      )
    } else {
      fit_n <- mgcv::gam(fml.n, data = long_data_s.n, family = "binomial")
      fit_d <- mgcv::gam(fml.d, data = long_data_s.d, family = "binomial")
    }
    #
    g_n_hazard <- sapply(seq_along(abar), function(j) {
      pX <- W.n; pX[,"bin_id"] <- j
      if (!is.null(SL.library)) as.numeric(predict(fit_n, pX)$pred)
      else as.numeric(predict(fit_n, type = "response", newdata = pX))
    })
    
    g_d_hazard <- sapply(seq_along(abar), function(j) {
      pX <- W.d; pX[,"bin_id"] <- j
      if (!is.null(SL.library)) as.numeric(predict(fit_d, pX)$pred)
      else as.numeric(predict(fit_d, type = "response", newdata = pX))
    })
    
    g_n_density <- t(apply(g_n_hazard, 1, hazard_to_density))
    g_d_density <- t(apply(g_d_hazard, 1, hazard_to_density))
    
    list(g_n = g_n_density, g_d = g_d_density)
  })
  #
  g.n <- lapply(result_list, `[[`, "g_n")
  g.d <- lapply(result_list, `[[`, "g_d")
  #
  return(list(g.n,g.d))
  #
}

as_ped<- function(formula, data, cut) {
  # Check inputs
  if (!inherits(formula, "formula")) stop("`formula` must be a formula.")
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (anyNA(cut) || !is.numeric(cut)) stop("`cut` must be a numeric vector without NAs.")
  #
  mf <- model.frame(formula, data)
  surv_obj <- mf[[1]]
  if (!inherits(surv_obj, "Surv")) stop("Left-hand side of formula must be a Surv object.")
  #
  y <- surv_obj[, "time"]
  status <- surv_obj[, "status"]
  XX <- mf[, -1, drop = FALSE]
  #
  intervals <- data.frame(
    start = head(c(-1, cut), -1),
    stop  = cut
  )
  #
  ped_list <- lapply(seq_along(y), function(i) {
    t_i <- y[i]
    X_i <- XX[i, , drop = FALSE]
    active_rows <- which(intervals$start < t_i)
    if (length(active_rows) == 0) return(NULL)
    
    int_i <- intervals[active_rows, , drop = FALSE]
    int_i$stop <- pmin(int_i$stop, t_i)
    
    # Add ped_status: 1 if event occurs in last interval
    int_i$ped_status <- 0
    if (status[i] == 1) {
      int_i$ped_status[nrow(int_i)] <- 1
    }
    
    cbind(
      id = i,
      bin_id = int_i$stop,
      ped_status = int_i$ped_status,
      X_i[rep(1, nrow(int_i)), , drop = FALSE]
    )
  })
  #
  ped_data <- do.call(rbind, ped_list)
  rownames(ped_data) <- NULL
  #
  return(ped_data)
}


hazard_to_density <- function(h) {
  d <- numeric(length(h))
  d[1] <- h[1]
  if (length(h) > 1) {
    for (k in 2:length(h)) {
      d[k] <- prod(1 - h[1:(k - 1)]) * h[k]
    }
  }
  return(d)
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

num_knots_generator <- function(max_degree, smoothness_orders, base_num_knots_0 = 500, 
          base_num_knots_1 = 200) 
{
  if (all(smoothness_orders > 0)) {
    return(sapply(seq_len(max_degree), function(d) {
      round(base_num_knots_1/2^(d - 1))
    }))
  }
  else {
    return(sapply(seq_len(max_degree), function(d) {
      round(base_num_knots_0/2^(d - 1))
    }))
  }
}

non.na.identical <- function(v1,v2){
  compare <- v1==v2
  all(na.omit(compare))
}

utils::globalVariables("i")

# not needed anymore, left for future
previous <-function(Qs,Ys,Ls,time){
  if(Qs[time]%in%Ys){
    Qlist <- Qs[Qs%in%Ys]
    pos <- match(Qs[time],Qlist)
    if(pos>1){prev<-Qlist[pos-1]}else{prev<-99}}else{
      Qlist <- Qs[Qs%in%Ls]
      pos <- match(Qs[time],Qlist)
      if(pos>1){prev<-Qlist[pos-1]}else{prev<-99}
  }
  return(prev)
}

