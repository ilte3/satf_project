get.bins <- function(x, 
                     rt_col = "rt[0-3][0-9]",
                     res_col = "key[0-3][0-9]",
                     window = c("fixed","free"),
                     offset = TRUE,
                     sauce = FALSE,
                     endbin = 16, binsize = 350, 
                     auditory = FALSE, audio.rt.base = 5006, 
                     quiet = TRUE,
                     plot = FALSE,
                     add.two = FALSE,
                     shift.edge = FALSE)
{
  if(!is.numeric(rt_col)) rt_col <- grep(rt_col,colnames(x), value=TRUE)
  if(!is.numeric(res_col)) res_col <- grep(res_col,colnames(x), value=TRUE)
  if(length(rt_col) != length(res_col)) stop("Number of RTs and Number of Responses did not match")
  ntrial <- nrow(x)
  bins <- (0:endbin) * binsize
  window <- match.arg(window)
  if(offset)
  {
    optbin <- opt.bins(na.omit(unlist(x[,rt_col])), endbin = endbin, type=window, add.two=add.two)
    if(!quiet) print(optbin)
    if(plot) plot(optbin) #plot if requested
    ofst <- optbin$offset
    bins <- bins + as.numeric(ofst) - 0.001 #not sure why -0.001 is needed, but is here anyway
  }
  #if(sauce) bins[0:2] <- 0:2 * binsize #no idea what sauce is for
  
  if(auditory){
    if(!("SentLength" %in% names(x))) stop("SentLength not found in x")
    if(!is.numeric(x$SentLength)) x$SentLength <- as.numeric(as.character(x$SentLength))
    x[,rt_col] <- x[,rt_col] - (x$SentLength - audio.rt.base) # this does the deriveRT thingy
  }
  
  ### need to flag if the two adjecent responses have the same RT
  flags <- data.frame(t(apply(x[,rt_col], 1, diff)))
  flags$rt31 <- NA
  flags[is.na(flags)] <- 1
  flags <- ifelse(flags == 0, 1, 0)
  colnames(flags) <- flag_col <- paste("flag",1:ncol(flags), sep=".")
  
  #create a data frame with required columns
  needed <-  c("Subject", "Session", "Block", "Set", "Condition", "Item", "CorrResp")
  wide <- cbind(x[,c(needed, rt_col, res_col)], flags)
  
  #chage it to long format, and remove NA
  long <- reshape(wide, varying=list(res_col, rt_col, flag_col), 
                  idvar="id", direction="long", v.names=c("RT","Response", "SeqFlag"))
  #long <- na.omit(long)
  colnames(long) <- c("Subject","Session","Trial","Experiment","Condition",
                      "Item","CorrResp","ResponseNumber","Response","Latency","SeqFlag", "id")
  rownames(long) <- NULL
  long$Accuracy <- as.numeric(long$Response == long$CorrResp)
  long$Bin <- findInterval(long$Latency, vec=bins)
  long$Bin[long$Bin %in% c(0,endbin+1)] <- 999;
  
  list(bins=long, opt.bins=optbin)
  #return(long)
}


opt.bins <- function(RT, maxRT=6300, breaksize=7, startbin=1,
                     endbin=16, binsize = 350, type=c("fixed","free"), add.two=FALSE, shift.edge = FALSE)
{
  ## first get the RT counts by breaksize and smooth
  rt.h <- hist(RT[RT<maxRT], breaks=seq(from=0, to=maxRT, by=breaksize), plot=FALSE)
  counts <- smooth(rt.h$counts,kind="3R")
  
  ## define hsw fitting function
  .fit.hsw <- function(midpoints, counts, amplitude, omega, phaseinterval){
    .hsw <-function(phase, midpoints, amplitude, omega){
      HSWvalue <- amplitude * sin((midpoints - phase) * (2*pi) / omega )
      HSWvalue[HSWvalue < 0] <- 0
      return(HSWvalue)
    }
    fitseries <- sapply(phaseinterval, FUN=.hsw, midpoints, amplitude, omega)
    fit <- apply(fitseries, 2, function(x,xp){sqrt(mean((x-xp)^2))}, counts)
    fit
  }
  
  ##
  binbreaks <- round(binsize/breaksize)
  if(shift.edge){
    firstbin <- startbin*binbreaks + (binbreaks + 1) # based on mrsatfnx.R
    lastbin <- endbin*binbreaks + (2 * binbreaks)   # based on mrsatfnx.R
  } else {
    firstbin <- startbin*binbreaks - binbreaks + 1   #based on Pyeongwhan's
    lastbin <- endbin*binbreaks                      #based on Pyeongwhan's
  }
  if(lastbin > length(counts)) stop("endbin too large for the given maxRT")
  bingroup <- factor(rep(startbin:endbin, each=binbreaks))
  peaks <- tapply(counts[firstbin:lastbin], bingroup, max)
  
  amp <- mean(peaks)
  
  type <- match.arg(type)
  if(type == "fixed")
  {
    omega <- ifelse(add.two, binsize + 2, binsize)
  } else {
    peaktime <- tapply(counts[firstbin:lastbin], bingroup, function(x){floor(median(which(x==max(x))))}) + 
      seq(firstbin,lastbin, by=binbreaks);
    omega <- median(diff(peaktime)) * breaksize 
  }
  
  by <- .5
  phaseinterval <- seq(from=0, to=(omega-by), by=by)
  d.fits <- .fit.hsw(rt.h$mids, rt.h$counts, amp, omega, phaseinterval);
  d.phase <- phaseinterval[which.min(d.fits)]
  d.offset <- abs(d.phase - omega/4);
  
  # now, create the res list
  res <- list(maxRT = maxRT,
              breaksize = breaksize,
              endbin = endbin,
              binsize = binsize,
              type = type, 
              offset = d.offset, 
              fit = min(d.fits),
              phase = d.phase,
              omega = omega,
              rt.hist=rt.h, 
              amp=amp
              )
  class(res) <- "optbin"
  res
}

print.optbin <- function(x, ...)
{
  cat(paste("Estimation type: ", x$type, sep=""))
  cat(paste("\nWindow size: ", round(x$omega), " ms", sep="")) 
  cat(paste('\nAmplitude estimate: ', round(x$amp, digits=2), sep=""));
  cat(paste("\nPhase: ",round(x$phase), sep=""));
  cat(paste("\nFit distance: ",round(x$fit, digits=2), sep=""));
  cat(paste("\nOffset: ",round(x$offset), "\n", sep=""));
}

plot.optbin <- function(x, 
                        main = "Response counts",
                        xlab="time (ms)", 
                        ylab = "response count",
                        col.hist = "dark grey", 
                        col.line = "red",
                        lty.line = "dotted",
                        pch = "+", 
                        cex=0.8,
                        legend = TRUE,
                        ...)
{
  rt.hist <- x$rt.hist
  plot(rt.hist, main=main, border=col.hist, xlab=xlab, ylab=ylab, xlim=c(0, x$maxRT), ...)
  abline(v=seq(from=x$offset, to=x$maxRT, by=x$omega), lty=lty.line, col=col.line)  
  
  .hsw <-function(x, midpoints, A, msOmega){
    HSWvalue <- A * sin((midpoints - x) * (2*pi) / msOmega )
    HSWvalue[HSWvalue < 0] <- 0;
    HSWvalue
  }
  points(rt.hist$mids, .hsw(x=x$phase, rt.hist$mids, x$amp,msOmega=x$omega), pch = pch, cex=cex)
  show.info <- ifelse(is.logical(legend), legend, TRUE)
  if(show.info)
  {
    if(is.logical(legend))
    {
      legend <- c(paste("Estimation type:", x$type),
                 paste("Window size:", round(x$omega), "ms"),
                 paste('Amplitude estimate:', round(x$amp)),
                 paste("Phase:",round(x$phase)),
                 paste("Fit distance:",round(x$fit, 3)),
                 paste("Offset:",round(x$offset)))
    }
    legend(0, y=max(rt.hist$count), legend=legend, cex=.7)
  }

}
