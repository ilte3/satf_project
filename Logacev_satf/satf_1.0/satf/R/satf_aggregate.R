
satf_aggregate_nyes <- function(data, id, time.id=c(), signal='signal', 
                                time='time', dv='response') 
{
  stopifnot.binary( data[[ dv[1] ]] )
  stopifnot(!(signal %in% id))
  stopifnot(!(signal %in% time.id))
  data <- ddply(data, c(id, signal), function(d) {
    d$n.responses.yes <- sum(d[[ dv[1] ]])
    d$n.responses <- length(d[[ dv[1] ]])
    for(colname in colnames(d)) {
      if(length(unique(d[[colname]])) > 1) {
        newcolname <- paste('mean', colname, sep='.')
        d[[ newcolname ]] <- mean( d[[colname]] )
        d[[colname]] <- NULL
      }
    }
    d[1,]
  })
  data[[ dv[1] ]] <- NULL
  mean.time <- paste('mean', time, sep='.')
  if(length(time.id) > 0 ) {
    data <- ddply(data, time.id, function(d) {
      if(time %in% colnames(d))
        d[[ time ]] <- mean(d[[ time ]])
      else
        d[[ mean.time ]] <- mean(d[[ mean.time ]])
      d
    })
  }
  data
}

satf_aggregate_dprime <- function(data, id, signal, dv=NULL, flat.min=0.5) {
  if(is.null(dv))
    dv <- c('n.responses.yes', 'n.responses')
  
  res <- ddply(data, id, function(d) {
    d.noise  <- d[ !as.logical(d[[signal]]), ]
    d.signal <- d[ as.logical(d[[signal]]), ]
    #reportifnot(nrow(d.signal) == 1, sprintf("ncol(d.signal)=%d", nrow(d.signal)))
    reportifnot(nrow(d.noise) == 1, sprintf("ncol(d.noise) = %d", nrow(d.noise)))
    
    hits = d.signal[[ dv[1] ]]
    misses = d.signal[[ dv[2] ]] - hits
    fas = d.noise[[ dv[1] ]]
    crs = d.noise[[ dv[2] ]] - fas
    d.signal[[ dv[1] ]] <- NULL
    d.signal[[ dv[2] ]] <- NULL
    
    summary <- ldply(1:nrow(d.signal), function(i) {
      compute_dprime(hits=hits[i], misses=misses[i], fas=fas, crs=crs, flat.min=flat.min)
    })
    cbind(d.signal, summary)
  })
  res
}



compute_dprime <- function(hits, fas, misses, crs, flat.min=.5)
{
  # 'correct' possible probabilites of 1  
  n.signal <- hits+misses; n.noise <- fas+crs
  if(n.signal < n.noise) {
    hits <- hits+flat.min; misses <- misses+flat.min;
    fas <- fas+flat.min*(n.noise/n.signal); crs <- crs+flat.min*(n.noise/n.signal);
  } else {
    hits <- hits+flat.min*(n.signal/n.noise); misses <- misses+flat.min*(n.noise/n.signal);
    fas <- fas+flat.min; crs <- crs+flat.min;
  }
  n.signal <- hits+misses; n.noise <- fas+crs
  
  p.hit <- hits/n.signal
  p.CR <- crs/n.noise
  dprime <- qnorm(p.hit)-qnorm(1-p.CR)
  c    <- -0.5*(qnorm(p.hit)+qnorm(1-p.CR))
  dprime.var <- ( p.hit*(1-p.hit) ) / ( n.signal*dnorm(qnorm(p.hit))^2 ) +
    ( p.CR*(1-p.CR) ) / ( n.noise*dnorm(qnorm(p.CR))^2 )
  
  c(dprime=dprime, dprime.var=dprime.var, c=c, c.var=dprime.var*0.25)
}


