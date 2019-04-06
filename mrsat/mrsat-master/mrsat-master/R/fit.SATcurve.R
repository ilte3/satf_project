fit.SATcurve <- function (data, par.cond, 
                          fix.asym = FALSE,
                          params = get.param(par.cond, auto.asym=fix.asym, data=data),
                          opt = c("acp" ,"nlminb", "nlm", "optim"), 
                          multi.opt = FALSE,
                          rep=1, quiet=TRUE, maxit=500, strict=TRUE, trace=FALSE)
{
  
  if(!is.data.frame(data) || !is.list(data)) data <- read.dat(data)
  
  npoints <- nrow(data)
  lags <- data$lags
  dprimes <- data$dprimes
  cond <- as.numeric(factor(data$condition))
  var.dp <- var(dprimes)
  
  if(missing(par.cond) || is.null(par.cond))
  {
    warning("par.cond is missing or null")
    par.cond <- list(asym=cond, rate=cond, incp=cond)
    params <- get.param(par.cond)
  }
  
  
  if(!is.data.frame(params)) params <- read.param(params)
  
  npar <- nrow(params)
  fixed.vals <- NULL
  if(fix.asym)
  {
    n.asym <- length(unique(par.cond$asym))
    fixed.vals <- params[1:n.asym,"start"]
    params <- params[(n.asym+1):npar,]
  }
   
  startvals <- params$start
  lower <- params$lower
  upper <- params$upper
  npar <- length(startvals)
  
  .optfun <- function(x, dprimes, lags, condition, par.cond=NULL, fixed.vals = NULL)
  {
    x <- c(fixed.vals, x)
    sum((dprimes - vareq(x, lags=lags, condition=condition, par.cond=par.cond))^2)
  }
  
  fit <- list()
  
  #options <- match.arg(opt, several.ok=T)
  options <- match.arg(opt, c("acp" ,"nlminb", "nlm", "optim"), several.ok = T)
  if(!multi.opt){ options <- options[1] }
  n.options <- length(options)
  
  if(rep > 2) options <- rep(options, times=rep)

  rep <- length(options)


  if(rep > 2 & !quiet) pb <- txtProgressBar(min = 1, max = rep, char = "=", style = 3)
  for(i in 1:rep)
  {
    if (rep > 2 & !quiet) setTxtProgressBar(pb, i)
    opt_method <- options[i]
    if(i > n.options) startvals <- runif(npar, min=lower, max=upper)
    
    if (opt_method == "acp")
    {
      optRes <- acp(start = startvals, objective = .optfun, 
                    dprimes=dprimes, lags=lags, condition=cond, par.cond=par.cond, fixed.vals=fixed.vals,
                    control = list(maxit=maxit, strict=strict), 
                    lower=lower, upper=upper
                    )
      optRes$estimate <- optRes$par
      optRes$par <- c(fixed.vals, optRes$par)
      optRes$method <- opt_method
    }
    else if (opt_method == "nlminb") 
    {
        optRes <- nlminb(start = startvals, objective = .optfun, 
                       lower=lower, upper=upper, 
                       control = list(iter.max=maxit, abs.tol=1e-20, rel.tol=1e-15), 
                       lags=lags, dprimes=dprimes, condition=cond, par.cond=par.cond, fixed.vals=fixed.vals)
        optRes$estimate <- optRes$par
        optRes$par <- c(fixed.vals, optRes$par)
        optRes$value <- optRes$objective
        optRes$method <- opt_method
        optRes$objective <- NULL
    }
    else if (opt_method == "nlm")
    {
        optRes <- nlm(p = startvals, f = .optfun,
                       typsize=lower, iterlim=maxit,
                       lags=lags, dprimes=dprimes, condition=cond, par.cond=par.cond, fixed.vals=fixed.vals)
        optRes$value <- optRes$minimum
        optRes$par <- c(fixed.vals, optRes$estimate)
        optRes$convergence <- optRes$code
        optRes$method <- opt_method
        optRes$minimum <- NULL
        
    }
    else if (opt_method == "optim"){
        optRes <- optim(par = startvals, fn = .optfun, method = "L-BFGS-B", 
                      lower=lower, upper=upper, 
                      control = list(maxit=maxit), 
                      lags=lags, dprimes=dprimes, condition=cond, par.cond=par.cond, fixed.vals=fixed.vals)
        optRes$estimate <- optRes$par
        optRes$par <- c(fixed.vals, optRes$par)
        optRes$method <- opt_method
    }
    
    fit[[i]] <- optRes
  }
  minfit <- which.min(sapply(fit, "[[", "value"))
  rep.out <- NULL
  if(rep > 1)
  {
    .par.print <- function(x)
    {
      res <- c(x$par, x$value)
      par.nam <- paste("par", 1:length(x$par), sep="")
      names(res) <- c(par.nam, "value")
      res
    }
    rep.out <- data.frame(t(sapply(fit, .par.print)))
    row.names(rep.out) <- make.names(sapply(fit, "[[", "method"), unique=TRUE)
  }
  bfit <- fit[[minfit]]
  if(!quiet) cat("\nBest Fit:", minfit, "opt=",bfit$method, "(convergence:",bfit$convergence,")","\n")
  fitted <- vareq(bfit$par,lags=lags,condition=cond, par.cond=par.cond)

  
  res <- list(fit = bfit, 
              SSE = SSE <- bfit$value,
              MSE = MSE <- SSE/npoints,
              RMSE = sqrt(MSE),
              #logLik = bfit$value/(2*var.dp) - log(sqrt(pi*var.dp)),
              logLik = (npoints + npoints * log(2 * pi) + npoints * log(SSE/npoints) ) / -2 , 
              npar = npar, 
              n = npoints, 
              R2 = R2 <- cor(dprimes, fitted)^2,
              adjR2 = 1 - (1-R2) *((npoints-1)/(npoints-npar-1)),
              data = data.frame(data, fitted=fitted),
              par.cond = par.cond,
              rep.out = rep.out
              )
  attr(res$logLik, "df") <- attr(res$MSE, "df") <- npar
  attr(res$logLik, "nobs") <- attr(res$MSE, "nobs") <- npoints
  class(res$logLik) <- class(res$MSE) <- "logLik"
  attr(res, "df") <- length(startvals)
  class(res) <- c("list","SATcurve")
  res
}

print.SATcurve <- function(x, digits = getOption("digits"), ...)
{
  cat('\nParameters:\n');
  print(x$fit$par, digits=digits);
  
  cat('\nSum of squared deviation:\n');
  print(x$fit$value, digits=digits)
  
  cat('\nR squared:\n');
  print(x$R2, digits=digits);
  
  cat('\nAdjusted R squared:\n');
  print(x$adjR2, digits=digits);
  
  cat('\nAIC:\n')
  print(AIC(x$logLik), digits=digits);
  
  cat('\nNumber of free parameters:\n')
  print(attr(logLik(x),"df"))
  
  cat('\nopt method:',  x$fit$method, '\n')
  cat('convergence:', x$fit$convergence, '\n')
  if(!is.null(x$fit$message))
    cat('message:', x$fit$message, '\n')

}

plot.SATcurve <- function(x, condition, xlim, ylim, npoints = 100, 
                          xlab = "Time", ylab="d\'",
                          legend=TRUE, labels = TRUE, lab.cex = 1,
                          residual = FALSE, pch, col, lwd = 1, ...)
{
  if (!inherits(x, "SATcurve")) stop("object not of class \"SATcurve\"")
  #get the data
  data <- x$data
  if(!is.factor(data$condition)) data$condition <- factor(data$condition)
  #set up parameters
  par.cond <- x$par.cond
  pl <- lapply(par.cond, unique)
  num.par <- sum(sapply(pl, length))
  if(length(x$fit$par) != num.par){ stop("length of unique elements of par.cond should be the same as length(par)")}
  pl <- relist(x$fit$par, skeleton=pl)
  
  par <- list( asym = pl$asym[par.cond$asym], 
               rate = abs(pl$rate[par.cond$rate]),
               incp = abs(pl$incp[par.cond$incp])
  )  
  par.tab <- do.call("cbind", par)
  row.names(par.tab) <- levels(data$condition)

  # get the subset
  if(!missing(condition))
  {
    data <- data[data$condition %in% condition,]
    data$condition <- factor(data$condition)
  }
  cond <- levels(data$condition)
  ncond <- length(levels(data$condition))
  
  #get plotting options
  if(missing(pch)) pch <- 1:ncond
  if(length(pch) < ncond) pch <- rep(pch, length.out=ncond)
  
  if(missing(col)) col <- 1:ncond
  if(length(col) < ncond) col <- rep(col, length.out=ncond)
  
  if(missing(xlim)){
    xlim <- range(data$lags)
    xlim <- c(floor(xlim[1]), ceiling(xlim[2]))
  }
  # whether to plot residuals or actual data points.
  if(residual)
  {
    if(missing(ylim)){ ylim <- c(-3, 3)}
    data$residual <- data$dprimes - data$fitted
    
    plot(x=data$lags, y=data$dprimes, type="n", xlab=xlab, ylab=ylab, 
        xlim=xlim, ylim=ylim, ...)
    abline(h = 0)
    d <- split(data, data$condition)
    for(i in 1:length(cond))
    {
      
      with(d[[i]], points(x=lags, y=residual, pch=pch[i], col=col[i]))
    }
  }else{
    if(missing(ylim)){ ylim <- c(-.5, 5)}

    
    plot(x=data$lags, y=data$dprimes, type="n", xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, ...)
    d <- split(data, data$condition)

    for(i in 1:length(cond))
    {
      with(d[[i]], points(x=lags, y=dprimes, pch=pch[i], col=col[i]))
      # get the values for the curve
      my.par <- par.tab[cond[i],]
      my.lags <- seq(from=min(xlim), to=max(xlim), length.out=npoints)
      my.cond <- rep(cond[i], npoints)
      my.value <- vareq(par=my.par, lags=my.lags, condition=my.cond, par.cond=list(asym=1, rate=1, incp=1))
      lines(x=my.lags, y=my.value, col=col[i], lwd=lwd)
    }
    
    if(legend)
    {
      
      par.lab <- apply(par.tab, 1, function(x){paste(round(x,2), collapse=", ")})[cond]
      par.lab <- paste(cond, " (", par.lab, ")", sep = "")
      my.leg <- legend(x=min(xlim), y=max(ylim), par.lab, col=col, lty=1, pch=pch)
    } else { my.leg <- list(rect=list(top=max(ylim), h=-.25) ) }
    show.info <- ifelse(is.logical(labels), labels, TRUE)
    if(show.info)
    {
      if(is.logical(labels))
      {
        labels <- list(paste("Adjusted R2:", round(x$adjR2, 4)),
                    paste('RMSE:', round(x$RMSE, 4))
        )
      }
      top <- my.leg$rect$top - my.leg$rect$h - .25
      text(x = min(xlim), y = seq(from=top, to=top - .25, length.out=2), pos = 4, labels, cex=lab.cex)
    }
    
  }

}

summary.SATcurve <- function(object, ...)
{
  if (!inherits(object, "SATcurve")) stop("object not of class \"SATcurve\"")
  par.cond <- object$par.cond
  par <- object$fit$par
  
  pl <- lapply(par.cond, unique)
  num.par <- sum(sapply(pl, length))
  if(length(par) != num.par){ stop("length of unique elements of par.cond should be the same as length(par)")}
  
  pl <- relist(par, skeleton=pl)
  
  asym <- pl$asym[par.cond$asym]
  rate <- abs(pl$rate[par.cond$rate])
  incp <- abs(pl$incp[par.cond$incp])
  
  res <- as.list(c(asym, rate, incp, object$R2, object$adjR2, object$SSE, object$RMSE, object$logLik, -2*logLik(object), AIC(object), BIC(object), object$npar, object$n  ))
  #, object$fit$method))
  names(res) <- c(paste("asym", 1:length(asym), sep=""), 
                  paste("rate", 1:length(rate), sep=""), 
                  paste("incp", 1:length(incp), sep=""), 
                  "R2","adjR2", "SSE", "RMSEfit", "logLik", "Deviance", "AIC", "BIC", "npar", "n")
  res$method <- object$fit$method
  as.data.frame(res)
}

SATsummary.list <- function(x)
{
  do.call("rbind", lapply(x, summary.SATcurve))
}

logLik.SATcurve <- function(object, useMSE = FALSE, ...)
{
  if(useMSE) res <- object$MSE
  else res <- object$logLik
  res
}

extractAIC.SATcurve <- function(fit, scale, k = 2, useMSE = FALSE, ...)
{
  loglik <- logLik.SATcurve(fit, useMSE=useMSE)
  edf <- attr(loglik, "df")
  c(edf, -2 * loglik + k * edf)
}

