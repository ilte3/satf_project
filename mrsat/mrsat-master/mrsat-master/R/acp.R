acp <- function(start, objective, ..., control = list(), lower, upper)
{
  
  obj <- quote(objective(a, ...))
  if(missing(lower)) stop('lower is missing')
  if(missing(upper)) stop('upper is missing')
  if (any(lower != -Inf) || any(upper != Inf)) {
    lower <- rep_len(as.double(lower), length(start))
    upper <- rep_len(as.double(upper), length(start))
  }

  alast <- a <- start
  vlster <- eval(obj)

  con <- list(trace = FALSE, maxit = 1000L, lmax = 20L, strict = FALSE, abs.tol = .Machine$double.eps, auto.correct = TRUE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  ## check to see if the start values and the upper/lower are not the same
  if(any(start == lower) || any(start == upper))
  {
    warnings("start contains the same value as lower/upper, which may lead to convergence difficulties.\n")
    if(con$auto.correct)
    {
      start[which(start == lower)] <- (lower[which(start == lower)] + upper[which(start == lower)] ) / 2
      start[which(start == upper)] <- (lower[which(start == upper)] + upper[which(start == upper)] ) / 2
    }      
  }
  
  
  astep <- abs(start/25)
  if(any(start == 0)) return(0)
  cond1 <- ((start + 2.1 * astep) > upper)
  astep[cond1] <- (upper[cond1] - start[cond1]) / 2.1
  cond2 <- ((start - 2.1 * astep) < lower)
  astep[cond2] <- (start[cond2] - lower[cond2]) / 2.1
  abase <- start - 3 * astep
  
  ## set up trace if it is TRUE
  if(con$trace)
  {
    trace <- matrix(NA, nrow=con$maxit * 2 * con$lmax * length(start), ncol=4+length(start)) 
    colnames(trace) <- c("istep", "ips", "l", paste("par[", 1:length(start), "]", sep=""), "value")
  } 

  counter = 1
  for(istep in 1:con$maxit)
  {
      for(ips in 1:2)
      {
        vlast <- vlster
        alast <- a
        for(l in 1:con$lmax)
        {
          for(i in 1:length(start))
          {
            vlster = 0
            #try five different values and get the best
            for(j in 1:5)
            {
              a[i] = abase[i] + j * astep[i]
              sersq <- eval(obj)
              if(j == 1 || sersq < vlster) 
              {
                vlster <- sersq
                jstore <- j
              }
            }
            a[i] = abase[i] + jstore * astep[i]
            
            
            cond1 <- (jstore-1)*(jstore-5)
            if(cond1 > 0) break
            
            astep[i] <- astep[i] / 3
            if(cond1 == 0)
              astep[i] <- astep[i] * 12
            
            
            if(a[i] + 2.1 * astep[i] > upper[i] )
              astep[i] <- (upper[i] - a[i]) / 2.1
            if(a[i] - 2.1 * astep[i] < lower[i] )
              astep[i] <- (a[i] - lower[i]) / 2.1
            
            abase[i] <- a[i] - 3 * astep[i]
            
            if(con$trace) trace[counter,] <- c(istep, ips, l, a, vlster)
            counter = counter + 1
          }
         
        }
        if(vlster >= vlast) break
      }
      
      sersq <- vlster
      astep <- a - alast
      aup <- upper - .1*astep
      alw <- lower + .1*astep
      
      alast <- a
      a <- a + astep
      a[a>aup] <- aup[a>aup]
      a[a<alw] <- alw[a<alw]
      vlster <- sersq
      sersq <- eval(obj) #vareq(par=a, dprime=dprime, lags=lags, condition=condition)
      
      a <- alast
      
      if(con$strict){ 
        fin.cond <- all(abs(astep) < .Machine$double.eps)
      } else {
        fin.cond <- any(abs(astep) == 0)
      }
      if(fin.cond){
        convergence = 0;
        break;
      }

      astep <- abs(astep)
      cond1 <- ((a + 2.1 * astep) > upper)
      astep[cond1] <- (upper[cond1] - a[cond1]) / 2.1
      cond2 <- ((a - 2.1 * astep) < lower)
      astep[cond2] <- (a[cond2] - lower[cond2]) / 2.1
      abase <- a - 3 * astep
  }
  
  if(istep == con$maxit) convergence = 1

  res <- list(par = a,
              convergence = convergence,
              iterations = istep,
              value = sersq,
              trace = ifelse(con$trace, as.data.frame(na.omit(trace)), con$trace)
              )
  res
}




