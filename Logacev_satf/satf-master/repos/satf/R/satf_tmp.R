
satf.group <- function(start, data, metric="logLik", multilevel=T, ...)
{
  stopifnot(metric %in% c("logLik", "logLikRaw"))
  library(plyr)
  init.params <- ddply(data, .(subject), function(d.cur) {
    print(paste("subject", d.cur$subject[1]))
    res <- satf(start=start, data=d.cur, metric=metric, lambda=~1, beta=~1, delta=~1)
    start <- res$par$satf.params
    res <- satf(start=start, data=d.cur, metric=metric, ...)
    res <- c(LL=res$value, unlist(res$par))
    res
  })
  if(!multilevel) {
    return(list(params=init.params, hyperparams=NULL))
  }
  init.params$LL.params <- NULL
  init <- init.params
  init$subject <- NULL
  init$LL <- NULL
  hyperparams <- data.frame()
  init.names <- colnames(init)
  print(init)
  
  for(i in 1:ncol(init)) {
    name <- colnames(init)[i] 
    short.name <- strsplit(name, split=".params.")[[1]][2]
    # exclude outlying values
    x <- remove.outliers(init[,name], factor=3)
    m <- mean(x); v <- var(x)/2
    sgn <- 0
    if(all(init[,name] <= 0)) sgn = -1
    else if(all(init[,name] >= 0)) sgn = 1
    if(m == 0) m = 0.1
    scale <- v/m; shape <- abs(m/scale)
    hyperparams <- rbind(hyperparams, data.frame(param=name, param.short=short.name,
                                                 hyperparam="scale", val=scale, constrained=sgn,
                                                 offset=0.001))
    hyperparams <- rbind(hyperparams, data.frame(param=name, param.short=short.name,
                                                 hyperparam="shape", val=shape+1, constrained=1,
                                                 offset=1))
    #hyperparams <- rbind(hyperparams, data.frame(param=name, param.short=short.name,
    #                                  hyperparam="m", val=median(x), constrained=sgn,
    #                                  offset=0))
    #hyperparams <- rbind(hyperparams, data.frame(param=name, param.short=short.name,
    #                                  hyperparam="sd", val=sd(x)/4, constrained=1,
    #                                  offset=0))
  }
  print(hyperparams)
  
  colnames(init.params) <- sapply(strsplit(colnames(init.params), split=".params."), function(x) x[length(x)])
  start.names <- (colnames(init.params))[-c(1:2)]
  start.params.bySubj <- dlply(init.params, .(subject),  function(d) unlist(d[1,-c(1:2)]) ) 
  
  transform.hyperparams <- function(hyperparams) within(hyperparams, {
    val[constrained!=0] <- val[constrained!=0] - offset[constrained!=0]
    val[constrained==0] <- val[constrained==0] - sign(val[constrained==0])*offset[constrained==0]
    val[constrained!=0] <- log(val[constrained!=0]*constrained[constrained!=0])
  })
  untransform.hyperparams <- function(hyperparams) within(hyperparams, {
    val[constrained!=0] <- constrained[constrained!=0]*exp(val[constrained!=0])
    val[constrained!=0] <- val[constrained!=0] + offset[constrained!=0]
    val[constrained==0] <- val[constrained==0] + sign(val[constrained==0])*offset[constrained==0]
  })
  
  Eval.SATF.MLV <- function(hyperparams.vec, params.indices=nrow(hyperparams.vec),
                            subj.params=F, plot.params=F) 
  {
    hyperparams$val[params.indices] <- hyperparams.vec
    hyperparams <- untransform.hyperparams(hyperparams)
    ###   
    #   print.params <- hyperparams$val
    #   names(print.params) <- hyperparams$param.short
    #   print("current params")
    #   print(print.params)
    ###    
    fn.hyperparams <- dlply(hyperparams, .(param), function(d) {
      if(all(d$hyperparam == c("scale","shape"))) {
        scale <- d$val[1]; shape <- d$val[2];
        fn <- function(x) dgamma(x, scale=scale, shape=shape, log=T)
      } else if(all(d$hyperparam == c("m","sd"))) {
        m <- d$val[1]; sd <- d$val[2];
        fn <- function(x) dnorm(x, mean=m, sd=sd, log=T)
      } else {
        stop("This distribution is not implemented.")
      }
    })
    hyperparams1 <- hyperparams$val[c(T,F)]
    hyperparams2 <- hyperparams$val[c(F,T)]
    parameter.sign <- sign(hyperparams1)
    fn.params.LL <- function(params=NULL, action="compute.LL") {
      if(action=="compute.LL") {
        params <- params$satf.params
        params <- params*sign(hyperparams1)
        hyperparams1 <- abs(hyperparams1)
        LLs <- dgamma(params, scale=hyperparams1, shape=hyperparams2, log=T)
        return( sum(LLs) )
      }
      else if(action=="compute.LLs") {
        params <- params$satf.params
        params <- params*sign(hyperparams1)
        hyperparams1 <- abs(hyperparams1)
        LLs <- dgamma(params, scale=hyperparams1, shape=hyperparams2, log=T)
        return( LLs )
      } else if(action=="return.hyperparams1") {
        return(hyperparams1)
      } else if(action=="return.hyperparams2") {
        return(hyperparams2)
      }
    }
    
    global.start <- start.params.bySubj[[ asc(d.cur$subject[1]) ]]
    fits <- ddply(data, .(subject), function(d.cur) {
      start <- start.params.bySubj[[ asc(d.cur$subject[1]) ]]
      #print("start")
      #print(start)
      #print("/start")
      
      wrong.sign <- which(parameter.sign != sign(start))
      start[wrong.sign] <- 0.01*parameter.sign[wrong.sign]
      ###startLL <- fn.params.LL(list(satf.params=cur.start),  action="compute.LLs")
      ###start.impossible <- which(abs(startLL)==Inf)
      ###start[start.impossible] <- start.impossible+rnorm(length(start.impossible), sd=5)
      ###print(start)
      res <- satf(start=start, data=d.cur, metric=metric,
                  fn.params.LL=fn.params.LL, optim.digits=NA, ...)
      if(is.null(res)) {
        stop("Failed inital evaluation even after resetting params.")
      }
      ###   print("param LL")
      ###   print( res$par$satf.params[2] )
      ###   print( fn.params.LL(res$par,  action="compute.LLs")[2] )
      return( c(LL=res$value, unlist(res$par)) )
      #print(res$par$satf.params-start)
    }, .parallel=T)
    
    #print("fits")
    #print(fits)
    if(nrow(fits) != nrow(init.params)) {
      return(-Inf)
    }
    if(plot.params) {
      print("params")
      print(sum(fits$LL.params))
      #print(round(fits$LL.params,2))
      p1 <- ggplot(data=fits, aes(x=satf.params.lambda, y=0))+geom_point()+ 
        stat_function(fun=function(x) exp((fn.hyperparams[['satf.params.lambda']])(x)))+
        scale_x_continuous(limits=c(0,10))
      p2 <- ggplot(data=fits, aes(x=satf.params.beta, y=0))+geom_point()+ 
        stat_function(fun=function(x) exp((fn.hyperparams[['satf.params.beta']])(x)) )+
        scale_x_continuous(limits=c(0,5))
      p3 <- ggplot(data=fits, aes(x=satf.params.delta, y=0))+geom_point()+ 
        stat_function(fun=function(x) exp((fn.hyperparams[['satf.params.delta']])(x)))+
        scale_x_continuous(limits=c(0,1))
      print(multiplot(p1,p2,p3))
    }
    
    #print( fits$LL )
    print( sum(fits$LL) )
    if(subj.params)
      res <- fits
    else
      res <- sum(fits$LL)
    return(res)
  }
  
  do.plot <- T  
  start.params <- function(hyperparams) {
    params <- Eval.SATF.MLV(hyperparams$val, subj.params=T, plot.params=do.plot)
    params <- dlply(params, .(subject),  function(d) {
      params <- unlist(d[1,-c(1:2,ncol(params))])
      names(params) <- start.names
      params
    })
    params
  }
  
  # determine new start parameters for participants
  hyperparams <- transform.hyperparams(hyperparams)
  indices.all <- seq(1, nrow(hyperparams), 1)
  start.params.bySubj <- start.params(hyperparams)
  
  # Fit parameters separately, first.
  for(i in seq(1, nrow(hyperparams), 2)) {
    indices <- c(i, i+1)
    print(paste("Fitting parameter ",(i+1)/2," of ",nrow(hyperparams)/2,".", sep=""))
    res <- optim(hyperparams$val[indices], Eval.SATF.MLV, params.indices=indices,
                 plot.params=F, control=list(fnscale=-1, maxit=100, reltol=1e-4))
    hyperparams$val[indices] <- res$par
    start.params.bySubj <- start.params(hyperparams)
    print(res)
  }
  
  print("Fitting all parameters, reltol=10e-4.")
  res <- optim(hyperparams$val[indices.all], Eval.SATF.MLV,
               params.indices=indices.all, 
               control=list(fnscale=-1, maxit=500, reltol=1e-4))
  print(res)
  hyperparams$val <- res$par
  start.params.bySubj <- start.params(hyperparams)
  convergence4 <- res$convergence 
  
  print("Fitting all parameters, reltol=10e-5.")
  res <- optim(hyperparams$val[indices.all], Eval.SATF.MLV,
               params.indices=indices.all, 
               control=list(fnscale=-1, maxit=500, reltol=1e-5))
  print(res)
  hyperparams$val <- res$par
  start.params.bySubj <- start.params(hyperparams)
  convergence5 <- res$convergence 
  
  print("Fitting all parameters, reltol=10e-6.")
  res <- optim(hyperparams$val[indices.all], Eval.SATF.MLV,
               params.indices=indices.all, 
               control=list(fnscale=-1, maxit=500, reltol=1e-6))
  print(res)
  hyperparams$val <- res$par
  start.params.bySubj <- start.params(hyperparams)
  convergence6 <- res$convergence 
  
  print("Fitting all parameters, reltol=10e-7.")
  res <- optim(hyperparams$val[indices.all], Eval.SATF.MLV,
               params.indices=indices.all, 
               control=list(fnscale=-1, maxit=500, reltol=1e-7))
  print(res)
  hyperparams$val <- res$par
  start.params.bySubj <- start.params(hyperparams)
  convergence7 <- res$convergence 
  
  print("Fitting all parameters, reltol=1e-8.")
  res <- optim(hyperparams$val[indices.all], Eval.SATF.MLV,
               params.indices=indices.all, 
               control=list(fnscale=-1, maxit=10^6, reltol=1e-8))
  print(res)
  hyperparams$val <- res$par
  convergence8 <- res$convergence 
  
  subj.params <- Eval.SATF.MLV(hyperparams$val, subj.params=T, plot.params=do.plot)
  hyperparams <- untransform.hyperparams(hyperparams)
  hyperparams$LL <- res$value  
  hyperparams$'convergence.1e-4' <- convergence4
  hyperparams$'convergence.1e-5' <- convergence5
  hyperparams$'convergence.1e-6' <- convergence6
  hyperparams$'convergence.1e-7' <- convergence7
  hyperparams$'convergence.1e-8' <- convergence8
  list(params=subj.params, hyperparams=hyperparams)
}



remove.outliers <- function(x, factor) {
  distance <- abs(mean(x)-x);
  x <- x[order(distance, decreasing=T)];
  distance <- sort(distance, decreasing=T);
  exclude <- c()
  for(i in 2:length(x)) {
    if(distance[i-1]/distance[i] > factor) {
      exclude <- c(exclude, i-1)
    } else {
      break;
    }
  }
  if(length(exclude) > length(x)) {
    stop("More than half the data were excluded as outliers.")
  }
  if(length(exclude)) {
    print(paste("Excluding:", paste(x[exclude], collapse="")))
    x <- x[(-1*exclude)]
  }
  x
}
