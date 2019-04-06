get.dprime <- function(x, signal.list, noise.list, is.index = FALSE, 
                       subset = (x$SeqFlag==0 & x$Experiment==1),
                       binmax=16, 
                       dprime.correction=c("extreme", "loglinear", "none"),
                       negative.dprime = FALSE,
                       bintail=0)
{
  #first, trim down the x to the set we need
  if(!missing(subset)) x <- subset(x, subset)
  x <- subset(x, x$Bin <= binmax)
  #check if the signal and noise list looks OK
  if(missing(signal.list)) stop("signal.list missing")
  if(missing(noise.list)) stop("noise.list missing")
  # moidfy the length(noise.list) if needed
  ncond <- length(signal.list)
  if(ncond != length(noise.list)) warning("length(sinal.list) and length(noise.list) are different ... adjusting noise.list")
  if(ncond < length(noise.list)) noise.list <- noise.list[1:ncond]
  if(ncond > length(noise.list)) noise.list <- rep(noise.list, length.out=ncond)
  # check the condition names
  nonames <- c(is.null(names(signal.list)), is.null(names(noise.list)))
  if(all(nonames)) names(signal.list) <- names(noise.list) <- 1:ncond
  if(nonames[2]) names(noise.list) <- names(signal.list)
  if(nonames[1]) names(signal.list) <- names(noise.list)
  if(!setequal(names(signal.list), names(noise.list))) stop("names(signal.list) and names(noise.list) are different")
  
  if(!is.index){
    sig.check <- setdiff(unlist(signal.list), unique(x$Condition))
    noise.check <- setdiff(unlist(noise.list), unique(x$Condition))
    if(length(sig.check) + length(noise.check) > 0)
    {
      not.found <- paste(sort(c(sig.check, noise.check)), collapse = "  ")
      stop("unable to find ", not.found, " in x$Condition. Please check the signal.list and/or noise.list")
    }
  }
  
  #get the dprime correction type
  dc <- match.arg(dprime.correction)
  #define averaging function
  if(dc == "extreme") {
    .average <- function(x){
      m <- mean(x)
      correct <- 1
      if(m == 1){
        m <- (length(x) - .5) / length(x)
        correct <- 2
      } 
      if(m == 0) {
        m <- .5 / length(x)
        correct <- 3
      }
      c(mean=m, correction=correct, length=length(x))
    }
  } else if(dc == "loglinear"){
    .average <- function(x){
      m <- (sum(x)+.5) / (length(x) + 1)
      c(mean=m, correction=4, length=length(x))
    }
  } else {
    .average <- function(x){
      c(mean=mean(x), correction=1, length=length(x))
    }
  }
  if(bintail > 0)
  {
    newmax <- binmax - bintail
    x[x$Bin>newmax, "Bin"] <- newmax
  }
  
  res <- list()
  for(i in names(signal.list))
  {
    sidx <- signal.list[[i]]
    nidx <- noise.list[[i]]
    if(is.index) sidx <- unique(sort(x$Condition))[sidx]
    if(is.index) nidx <- unique(sort(x$Condition))[nidx]
    #return(x)
    hit <- aggregate(Accuracy ~ Bin, data=x, subset=x$Condition %in% sidx, .average)
    hit <-  data.frame(cbind(Bin=hit[,1], hit[,2]), stringsAsFactors = FALSE)
    #return(hit)
    fa <- aggregate((1 - Accuracy) ~ Bin, data=x, subset=x$Condition %in% nidx, .average)
    fa <- data.frame(cbind(Bin=as.numeric(fa[,1]), fa[,2]), stringsAsFactors = FALSE)

    tmp <- merge(hit, fa, by=c("Bin"), all=TRUE)

    names(tmp) <- c("bin","hit", "hit.correction", "hit.denom", "fa", "fa.correction", "fa.denom")
    tmp$lags <- aggregate( Latency/1000 ~ Bin, data=x, subset=x$Condition %in% c(sidx, nidx), mean)$`Latency/1000`

    tmp$dprimes <- qnorm(tmp$hit) - qnorm(tmp$fa)
    tmp$condition <- i
    res[[i]] <- tmp
  }
  res <- as.data.frame(do.call(rbind.data.frame, res))
  
  # negative dprime
  if(!negative.dprime) res$dprimes[res$dprimes < 0] <- 0

  errors <- c("none", "extreme(1)", "extreme(0)", "loglinear")
  res$hit.correction <- errors[res$hit.correction]
  res$fa.correction <- errors[res$fa.correction]
  #res$bin
  row.names(res) <- NULL
  if(any(is.infinite(res$dprimes))) warning("dprimes contain Inf/-Inf.")
  attr(res, "dprime.correction") <- dprime.correction
  class(res) <- c("mrsat.data", "data.frame")
  res
}

average.dprime <- function(bin.list)
{
  all <- do.call(rbind, bin.list)
  mean <- aggregate(cbind(lags, dprimes, hit, fa) ~ bin + condition, data=all, mean)
  se <- aggregate(cbind(lags, dprimes, hit, fa) ~ bin + condition, data=all, function(x){sd(x)/sqrt(length(x))})
  colnames(se)[3:6] <- c("lags.se", "dprimes.se", "hit.se", "fa.se")
  res <- merge(mean, se, by=c("condition", "bin"), sort=FALSE)
  class(res) <- c("mrsat.data","data.frame")
  res
}


plot.mrsat.data <- function(x, condition, xlim, ylim, legend=TRUE, pch, col, ...)
{
  if(!inherits(x, "mrsat.data")) stop("x not of class \"mrsat.data\"")
  data <- x
  if(!is.factor(data$condition)) data$condition <- factor(data$condition, levels =unique(data$condition), ordered = T)
  
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
  if(missing(ylim)){ ylim <- c(-.5, 5)}
  
  plot(x=data$lags, y=data$dprimes, type="n", xlab="Time", ylab="d\'", ylim=ylim, xlim=xlim, ...)
  d <- split(data, data$condition)
  
  for(i in 1:length(cond))
  {
    with(d[[i]], points(x=lags, y=dprimes, pch=pch[i], col=col[i]))
    with(d[[i]], lines(x=lags, y=dprimes, col=col[i]))
  }
  if(legend) legend(x=min(xlim), y=max(ylim), cond, col=col, pch=pch)
  
}
