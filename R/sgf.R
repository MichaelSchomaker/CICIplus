# Task: simulation 13 for survival with weights
# pass on options for calc.weights: i.e. parallelize for hal_density and screening options for model.formulas.update
# calc.weights: abar = matrix not possible so far; either stop() or make it possible
# survival: deterministic only when Q=Y, i.e. when Q is not a L-node
# individual interventions: if two schemes have identical "names", intervention assignment is not clear -> "check function"; Han
# categorical interventions? -> or categorical simply numeric and document?
# check function for correct survival setup?

sgf <- function(X, Anodes, Ynodes, Lnodes = NULL, Cnodes = NULL,
                abar =  NULL, survivalY = FALSE,
                SL.library = "SL.glm", SL.export = NULL,
                Yweights = NULL,
                calc.support = FALSE, B = 0, 
                ncores=1, verbose=TRUE, seed=NULL, prog=NULL,...){
  
  ### checks and setup ###
  if(is.null(seed)==FALSE){set.seed(seed)}
  if(is.null(prog)==FALSE){write(matrix("started with setup..."),file=paste(prog,"/progress.txt",sep=""))}
  model.families <- CICI:::assign.family(X)
  catint <- FALSE
  #
  if(is.data.frame(X)==FALSE){stop("'X' (i.e. your data) needs to be a data.frame")}
  if(missing(Anodes)){stop("'Anodes' is missing. Please specify your intervention node(s).")}  
  if(missing(Ynodes)){stop("'Ynodes' is missing. Please specify your outcome node(s).")}
  if(min(which(colnames(X)%in%Ynodes))<min(which(colnames(X)%in%Anodes))){stop("Ynodes can not occur before Anodes.\n Likely you specified pre-intervention variables as Ynode?")}
  if(is.character(Anodes)==FALSE & is.null(Anodes)==FALSE){stop("'Anodes' needs to be a character vector containing the respective variable names.")}
  if(is.character(Cnodes)==FALSE & is.null(Cnodes)==FALSE){stop("'Cnodes' needs to be a character vector containing the respective variable names.")}
  if(is.character(Lnodes)==FALSE & is.null(Lnodes)==FALSE){stop("'Lnodes' needs to be a character vector containing the respective variable names.")}
  if(is.character(Ynodes)==FALSE & is.null(Ynodes)==FALSE){stop("'Ynodes' needs to be a character vector containing the respective variable names.")}
  if(any(substr(colnames(X),1,4)=="list")){stop("Variable names that start with 'list' are not allowed. Please rename.")}
  if(any(c(Ynodes,Lnodes,Anodes,Cnodes)%in%colnames(X)==FALSE)){stop(paste("You have specified the following variable name(s) which are not part of the data:",
                                                                           paste(c(Ynodes,Lnodes,Anodes,Cnodes)[c(Ynodes,Lnodes,Anodes,Cnodes)%in%colnames(X)==FALSE],collapse=" ") ,"\n"))}
  if(max(which(colnames(X)%in%Anodes))==ncol(X)){stop("Anodes/Cnodes should not be in last column")}
  if(any(sapply(X,is.ordered))){stop("Ordered variables are not allowed")}
  if(is.logical(survivalY)==FALSE){stop("survivalY needs to be TRUE or FALSE")}
  if(survivalY==TRUE){if(verbose==TRUE){if(is.null(Cnodes)){cat("Warning: you indicate that you have survival data (survivalY=T), but you have no Cnodes specified. \n")}}}
  if(is.logical(verbose)==FALSE){stop("'verbose' needs to either TRUE or FALSE")}
  if(is.numeric(ncores)==FALSE){stop("'ncores' needs to be numeric")}
  if(is.character(prog)==FALSE & is.null(prog)==FALSE){stop("'prog' needs to be a character vector")}
  if(any(abar=="natural")){stop("Natural course scenario not possible for sequential g-formula. Use gformula().")}
  #
  if(any(model.families=="binomial")){bin.problem <- !sapply(subset(X,select=(model.families=="binomial")), CICI:::is.binary)
  if(any(bin.problem)){if(verbose==TRUE){cat(paste("Binary variables have been recoded:",paste(names(bin.problem)[bin.problem],collapse=","),"\n"))}
    X[,names(bin.problem)[bin.problem]] <- data.frame(sapply(subset(X,select=names(bin.problem)[bin.problem]),CICI:::binary.to.zeroone,verb=verbose)) }
  }
  #
  # find Q regression information
  loop.Q    <- find.Qs(dat=X, L=Lnodes, Y=Ynodes, A=Anodes, C=Cnodes) 
  Qnodes    <- unique(unlist(loop.Q[[1]]))
  model.Q.families <- model.families[colnames(X)%in%c(Qnodes)]
  if(any(substr(model.Q.families,1,4)=="mult")){stop("Categorical outcomes (or categorical Lnodes that act as outcome in Q-models) with >2 categories currently not supported")}
  if(any(model.Q.families=="poisson")){model.Q.families[model.Q.families=="poisson"]<-"gaussian"
  if(verbose==TRUE){cat("Note: your outcome (or Lnode that acts as outcome in Q-regression) is treated as continuous/gaussian in the Super Learner.\n")}}
  if(is.null(Yweights)==FALSE & any(model.Q.families=="binomial")){alert<-TRUE}else{alert<-FALSE} 
  if(is.null(Yweights)==TRUE  & any(model.Q.families=="binomial") & verbose==TRUE){cat("Note: with binomial outcomes (and >2 time points), the conditional outcome models need to model proportional data. \n Make sure your learners can deal with this.\n")}
  if(is.null(Yweights)==FALSE){if(is.list(Yweights)==FALSE | length(Yweights)!=max(unlist(lapply(loop.Q$Qs.for.each.Y,length))) | any(dim(Yweights[[1]])!=c(nrow(X),length(abar))) ){stop("Yweights do not have right format. Type ?sgf and read the details section.")}}
  
  ### matrices to store results ###
  n.a <- length(Anodes); n.t <- length(Ynodes); n.l <-length(Lnodes); time.points <- 1:n.t
  if(is.matrix(abar)==TRUE){interventions <- do.call(rbind, replicate(n.t, abar, simplify=FALSE)); i.type <- "custom"}else{
    if(is.vector(abar)){interventions <- do.call(rbind, replicate(n.t, matrix(rep(abar,n.a),ncol=n.a), simplify=FALSE)); i.type<-"standard"}
    if(is.list(abar)){interventions <- do.call(rbind,replicate(n.t,do.call(rbind,lapply(lapply(abar,colnames),as.numeric)),simplify=FALSE)) ; i.type="individual"}                          
    }
  n.int <- dim(interventions)[1]/n.t
  if(i.type=="custom" & calc.support==TRUE){calc.support<-FALSE; if(verbose==TRUE){cat("Note: no support can be calculated for custom intervention strategies.\n")}}
  if(i.type=="individual" & calc.support==TRUE){calc.support<-FALSE; if(verbose==TRUE){cat("Note: no support can be calculated for individual intervention strategies.\n")}}
  
  store.results           <- as.data.frame(matrix(NA,nrow=(n.t)*n.int,ncol=1+n.a+1))
  colnames(store.results) <- c("time",paste("a",1:n.a ,sep=""),"psi")
  store.results$time      <- rep(time.points ,each=n.int)
  store.results[,2:(n.a+1)] <-  interventions
  
  if(length(Ynodes)==length(Anodes)){needed <- !(store.results$time<matrix(rep(time.points,nrow(store.results)),nrow=nrow(store.results),byrow=T))}else{
    needed <-  matrix(which(colnames(X)%in%Anodes),nrow=n.int*length(time.points),ncol=length(which(colnames(X)%in%Anodes)),byrow=T) <
      matrix(rep(which(colnames(X)%in%Ynodes),each=n.int),nrow=n.int*length(time.points),ncol=n.a)  
  }
  store.results[,2:(n.a+1)][needed==FALSE] <- NA
  #
  SL.summary <- matrix(NA,ncol=length(Qnodes),nrow=length(SL.library),dimnames = list(unlist(lapply(SL.library,paste,collapse="_")),Qnodes))
  #
  
  ### Parallelization & Setup
  if(ncores>1){
    if(ncores > parallel::detectCores()){
      ncores <- parallel::detectCores();if(verbose==TRUE){cat(paste("Note: You only have",ncores,"threads which can be utilized. \n"))}
    }
    if(verbose==TRUE){cat(paste("Note: You initialized parallel computation using",ncores,"threads...initializing cluster now... \n"))}
    cl <- parallel::makeCluster(ncores); doParallel::registerDoParallel(cl)
    exp.var <- c(unlist(SL.library),SL.export,"multiResultClass","non.na.identical")
  }else{exp.var=NULL; foreach::registerDoSEQ()}
  
  if(B>0){analysis.b<- rep(list(NA),B)}else{analysis.b<-NULL}
  if(is.null(prog)==FALSE){if(verbose==TRUE){cat(paste("Progress will be saved in:",prog,"\n"))}
    write(matrix("started with sequential g-formula calculations in original data...\n"),file=paste(prog,"/progress.txt",sep=""),append=TRUE)}
  
  # Prepare support measures, if relevant
  if(calc.support==TRUE){
    support <- CICI:::calculate.support(dat=X,A=Anodes,intervention=interventions[store.results$time==1,])
    updat.index2 <- unlist(lapply(apply(t(outer(which(colnames(X)%in%Anodes), which(colnames(X)%in%Ynodes), "<")),1,which),max))   
  }
  
  ### ANALYSIS ###
  
  analysis <- foreach(i = 1:nrow(store.results), .export=exp.var, .packages=c("SuperLearner"),
                      .errorhandling="pass") %dorng% try({  
                        #
                        pind <- grep("predict",exp.var)
                        if(length(pind)>0){pf <- exp.var[pind]
                        for(j in 1:length(pf)){.S3method("predict", strsplit(pf[j],split="predict.")[[1]][2], get(pf[j]))}
                        }
                        #                
                        mydat <- X
                        current.t <- store.results[i,1]
                        current.Y <- Ynodes[current.t]
                        mydat <- mydat[,1:(which(colnames(mydat)==current.Y))]
                        ind <- i - (nrow(store.results)/n.t*(store.results[i,1]-1))
                        Q.info <- loop.Q[[1]][[current.t]]
                        
                        # Step 1: intervene
                        gdata <- mydat
                        all.rel.Anodes  <-   which(colnames(mydat)%in%c(Anodes))
                        if(i.type!="individual"){replacement <- matrix(rep((store.results[i,2:(n.a+1)])[is.na(store.results[i,2:(n.a+1)])==F], nrow(mydat)), nrow=nrow(mydat), byrow=T)}else{
                                                 replacement <- abar[[which(sapply(lapply(lapply(abar, colnames), as.numeric),
                                                                                   non.na.identical, as.vector(unname(unlist(store.results[i,2:(n.a+1)]))) ) ) ]]  
                                                 replacement <- replacement[,1:length(all.rel.Anodes)]
                                                 }
                        if(length(all.rel.Anodes)==1){replacement<-as.vector(replacement)}; 
                        gdata[,all.rel.Anodes] <-    replacement

                        # Step 2: sequential g-formula
                        for(t in length(Q.info):1)try({
                          if(is.null(Yweights)==FALSE){w.t<-Yweights[[t]][,ind]}else{w.t<-rep(1,nrow(X))} 
                          if(t==length(Q.info)){Y.t<-mydat[,which(colnames(mydat)==Q.info[length(Q.info)])]}else{Y.t <- Q.t}
                          if(t>1){Y.w <- Y.t*w.t}else{Y.w<-Y.t}
                          fdata <- cbind(mydat[,1:((which(colnames(mydat)==loop.Q[[2]][[current.t]][t])))],Y.w)
                          selind<-!is.na(fdata$Y.w);fdata <- fdata[selind,]
                          if(is.null(Cnodes)==FALSE){if(survivalY==TRUE){
                              incl <- setdiff(names(fdata),colnames(fdata)[colnames(fdata)%in%c(Cnodes,Ynodes)])}else{
                              incl <- setdiff(names(fdata),colnames(fdata)[colnames(fdata)%in%c(Cnodes)])}
                              incl2<- incl[incl!="Y.w"];fdata<-fdata[, incl]}else{incl2<-colnames(fdata)[colnames(fdata)!="Y.w"]}
                          if(alert==FALSE){suitable.family<-model.Q.families[Qnodes%in%Q.info[t]]}else{if(all(fdata$Y.w%in%c(0,1))){suitable.family<-"binomial"}else{suitable.family<-"gaussian"}}
                          m.Y   <- try(SuperLearner::SuperLearner(Y=fdata$Y.w, X=fdata[,-grep("Y.w",colnames(fdata))],
                                                                  SL.library=SL.library, family=suitable.family, ...), silent=TRUE) 
                          SL.summary[,colnames(SL.summary)%in%Q.info[t]] <- m.Y$coef
                          if(t>1){Q.t <- mydat[,Q.info[t-1]]}else{Q.t<- mydat[,Q.info[t]]} # implies: Q_t=1 if Y_t-1 = 1 
                          if(is.null(Cnodes)==FALSE & survivalY==FALSE){selind <- !is.na(Q.t)}
                          Q.t[selind] <- try(predict(m.Y, newdata = gdata[selind,incl2])$pred, silent=TRUE)
                          if(t==1){results<-weighted.mean(Q.t,w=w.t)}
                          #
                          if(model.Q.families[t]=="binomial" & is.null(Yweights)==TRUE){if(any(m.Y$library.predict<0)){
                            mes <- paste("Caution: negative predictions in:",paste(colnames(m.Y$library.predict)[apply(apply(m.Y$library.predict,2,function(x) x<0),2,any)],collapse = " ")); if(verbose==TRUE){print(mes)}
                            if(is.null(prog)==FALSE){try(write(matrix(mes),file=paste(prog,"/progress.txt",sep=""),append=T))}  }}
                          #
                        },silent=TRUE)
                        if(class(m.Y)=="try-error"){results<-NA}
                        #
                        all.results <- multiResultClass(); all.results$result1 <- results; all.results$result2 <- SL.summary
                        return(all.results)
                      },silent=TRUE)
  
  ##########
  # Step 5: Bootstrapping
  if(B>0){
    if(verbose==TRUE){cat("starting with bootstrapping \n")};if(is.null(prog)==FALSE){write(matrix("started with bootstrapping...\n"),file=paste(prog,"/progress.txt",sep=""),append=TRUE)}
    if(is.null(seed)==FALSE){set.seed(seed)}
    rng <- rngtools::RNGseq(B*nrow(store.results), seed); r <- NULL
    b.index <- apply(matrix(rep(1:nrow(X),B),ncol=B), 2, sample, replace=TRUE)
    boots <- foreach(b = 1:B) %:%
      foreach(i = 1:nrow(store.results), r=rng[(b-1)*nrow(store.results) + 1:nrow(store.results)],
              .export=exp.var, .packages=c("SuperLearner"), .errorhandling="pass") %dopar% {
                if(is.null(seed)==FALSE){rngtools::setRNG(r)}
                #
                pind <- grep("predict",exp.var)
                if(length(pind)>0){pf <- exp.var[pind]
                for(j in 1:length(pf)){.S3method("predict", strsplit(pf[j],split="predict.")[[1]][2], get(pf[j]))}
                }
                #   
                mydat <- X[b.index[,b],]
                if(is.null(prog)==FALSE & i==1){try(write(matrix(paste("performing calculations on bootstrap sample",b)),file=paste(prog,"/progress.txt",sep=""),append=T))}
                if(verbose==TRUE){if(i==1){cat(paste0("...",b));if(b%%10==0){cat("\n")} }}
                #
                current.t <- store.results[i,1]
                current.Y <- Ynodes[current.t]
                mydat <- mydat[,1:(which(colnames(mydat)==current.Y))]
                ind <- i - (nrow(store.results)/n.t*(store.results[i,1]-1))
                Q.info <- loop.Q[[1]][[current.t]]
                # Step 1: intervene
                gdata <- mydat
                all.rel.Anodes  <-   which(colnames(mydat)%in%c(Anodes))
                if(i.type!="individual"){replacement <- matrix(rep((store.results[i,2:(n.a+1)])[is.na(store.results[i,2:(n.a+1)])==F], nrow(mydat)), nrow=nrow(mydat), byrow=T)}else{
                  replacement <- abar[[which(sapply(lapply(lapply(abar, colnames), as.numeric), non.na.identical, as.vector(unname(unlist(store.results[i,2:(n.a+1)]))) ) ) ]]  
                  replacement <- replacement[,1:length(all.rel.Anodes)]
                }
                if(length(all.rel.Anodes)==1){replacement<-as.vector(replacement)}; 
                gdata[,all.rel.Anodes] <-    replacement
                
                # Step 2: sequential g-formula
                for(t in length(Q.info):1)try({
                  if(is.null(Yweights)==FALSE){w.t<-Yweights[[t]][,ind]}else{w.t<-rep(1,nrow(X))} 
                  if(t==length(Q.info)){Y.t<-mydat[,which(colnames(mydat)==Q.info[length(Q.info)])]}else{Y.t <- Q.t}
                  if(t>1){Y.w <- Y.t*w.t}else{Y.w<-Y.t}
                  fdata <- cbind(mydat[,1:((which(colnames(mydat)==loop.Q[[2]][[current.t]][t])))],Y.w)
                  selind<-!is.na(fdata$Y.w);fdata <- fdata[selind,]
                  if(is.null(Cnodes)==FALSE){if(survivalY==TRUE){
                    incl <- setdiff(names(fdata),colnames(fdata)[colnames(fdata)%in%c(Cnodes,Ynodes)])}else{
                    incl <- setdiff(names(fdata),colnames(fdata)[colnames(fdata)%in%c(Cnodes)])}
                    incl2<- incl[incl!="Y.w"];fdata<-fdata[, incl]}else{incl2<-colnames(fdata)[colnames(fdata)!="Y.w"]}
                  if(alert==FALSE){suitable.family<-model.Q.families[Qnodes%in%Q.info[t]]}else{if(all(fdata$Y.w%in%c(0,1))){suitable.family<-"binomial"}else{suitable.family<-"gaussian"}}
                  m.Y   <- try(SuperLearner::SuperLearner(Y=fdata$Y.w, X=fdata[,-grep("Y.w",colnames(fdata))],
                                                          SL.library=SL.library, family=suitable.family, ...), silent=TRUE) 
                  if(t>1){Q.t <- mydat[,Q.info[t-1]]}else{Q.t<- mydat[,Q.info[t]]} # implies: Q_t=1 if Y_t-1 = 1 
                  if(is.null(Cnodes)==FALSE & survivalY==FALSE){selind <- !is.na(Q.t)}
                  Q.t[selind] <- try(predict(m.Y, newdata = gdata[selind,incl2])$pred, silent=TRUE)
                  if(t==1){results<-weighted.mean(Q.t,w=w.t)}
                },silent=TRUE)
                if(class(m.Y)=="try-error"){results<-NA}
                #############################################
                return(results)
              }
    
  }
  ##########
  
  if(ncores>1){parallel::stopCluster(cl)}
  
  if(length(SL.library)>1){SL.summary <- apply(simplify2array(lapply(analysis, '[[', 2)), 1:2, mean, na.rm=T)}else{
    SL.summary <- matrix(1); rownames(SL.summary)<-SL.library
  }
  analysis   <- do.call("rbind",lapply(analysis, '[[', 1))
  store.results$psi <- c(analysis)
  
  if(B>0){
    for(b in 1:B){analysis.b[[b]] <- do.call("rbind",lapply(boots[[b]], '[[', 1))}
    boot.failure <- lapply(analysis.b,is.character)
    if(sum(unlist(boot.failure))>0){boots <- boots[-c(1:B)[unlist(boot.failure)]];analysis.b <- analysis.b[-c(1:B)[unlist(boot.failure)]]
    if(verbose==TRUE){cat(paste("Caution:",sum(unlist(boot.failure)),"bootstrap sample(s) were removed due to errors \n"))}   }
    newB <- B-sum(unlist(boot.failure))
    store.results[,c("l95","u95")] <-  t(apply(matrix(unlist(analysis.b),ncol=newB),1,quantile,probs=c(0.025,0.975)))
  }
  
  
  # calculate support if desired
  if(calc.support==TRUE){
    if(n.int<6 & n.int>2 & verbose==TRUE){cat(paste("Note: you have specified only",n.int,"interventions. The reported support diagnostics may not be reliable here. \n"))}
    diagn <- list(crude_support=support$crude_support,conditional_support=support$cond_support)
    cn <- paste("a",1:n.a ,sep=""); if(i.type=="standard"){rn <- as.character(abar)}else{rn <- paste("Strategy",1:nrow(abar))}
    diagn <- lapply(diagn,as.data.frame)
    rownames(diagn$crude_support) <- rn; rownames(diagn$conditional_support) <- rn; colnames(diagn$crude_support) <- cn; colnames(diagn$conditional_support) <- cn
  }else{diagn<-NULL}
  
  if(is.null(prog)==FALSE){write(matrix("finished calculations.\n"),file=paste(prog,"/progress.txt",sep=""),append=TRUE)}
  
  # return results
  res= list(results=store.results, 
            diagnostics=diagn, 
            SL.weights = round(SL.summary, digits=2),
            setup=list(i.type = i.type, n.t=n.t, B=B, fams=model.families, measure="default",
                       Ynodes = Ynodes, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, abar=abar,
                       support=calc.support, survival=survivalY, Qblocks=loop.Q, Qnodes=Qnodes, catint=catint)
  )
  
  class(res) <- "gformula"
  res
  
  #
}