vareq <- function(par, lags, condition, 
                  par.cond=list(asym=condition, rate=condition, incp=condition))
{
  
  if(!is.numeric(condition))
    condition <- as.numeric(as.factor(condition))
  
  cond <- unique(condition)
  n.cond <- length(cond)
  
  if (is.matrix(par.cond)){
      pdim <- dim(par.cond)
      test <- pdim != c(length(cond), 3)
      if(any(test)) stop(paste("dimension of par.cond (", pdim[1], ",", pdim[2], ") is incorrect!", sep=""))
      par.cond <- list(
        asym = c(par.cond[,1])[condition],
        rate = c(par.cond[,2])[condition],
        incp = c(par.cond[,3])[condition]
      )
  }

  if(is.list(par.cond))
  {
    if(length(par.cond) != 3) stop("length(par.cond) must be 3")
    if(all(lapply(par.cond, length) == n.cond))
    {
      par.cond <- list(
        asym = par.cond$asym[condition],
        rate = par.cond$rate[condition],
        incp = par.cond$incp[condition]
        )
    }
    else if(all(lapply(par.cond, length) != length(condition)))
    {
      stop("length of each element of par.cond should be the same as either the number of conditions or length(condition)")
    }
    
    
    pl <- lapply(par.cond, unique)
    num.par <- sum(sapply(pl, length))
    if(length(par) != num.par){ stop("length of unique elements of par.cond should be the same as length(par)")}
    
    pl <- relist(par, skeleton=pl)

  }

  asym <- pl$asym[par.cond$asym]
  rate <- abs(pl$rate[par.cond$rate])
  incp <- abs(pl$incp[par.cond$incp])
  vareq <- asym * (1.0 - exp( -rate * (lags - incp) ))
  vareq[lags < incp] <- 0
  vareq
}
