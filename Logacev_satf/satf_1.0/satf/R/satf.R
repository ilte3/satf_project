.packageName <- "satf"

Deviance <- function(LL) -2*LL
AIC <- function(LL, k) Deviance(LL) + 2*k
BIC <- function(LL, k, n) Deviance(LL) + k*log(n)

p2logodds <- function(p) log(p/(1-p))
logodds2p <- function(lodds) exp(lodds)/(1+exp(lodds))

NumSmallest <- -.Machine$double.xmax
NumLargest <- .Machine$double.xmax


compute_logLikFn <- function(coefs, by_row=FALSE, force_update=FALSE, invert=FALSE) {
	res = rcpp_compute_logLikFn(coefs, by_row, FALSE, force_update)
#  print(res)
  if(invert) return(-res)
  else       return(res)
}

compute_logLikFn_gradient <- function(coefs, by_row=FALSE) {
  res = rcpp_compute_logLikFn_gradient(coefs, by_row, FALSE)
  colnames(res) = names(coefs)
#  print(rbind(coefs,res))
  res
}

compute_logLikFn_gradient_BHHH <- function(coefs) {
  res = rcpp_compute_logLikFn_gradient(coefs, TRUE, FALSE)
  colnames(res) = names(coefs)
#  print(rbind(coefs,res))
#  stop()
  res
}

parameter_combinations <- function(namedlist) {
  membernames <- names(namedlist)
  cross_prod = cset(namedlist[[1]])
  for(i in 2:length(namedlist))
    cross_prod = cross_prod*cset(namedlist[[i]])
  llply(as.list(cross_prod), function(x) { x= unlist(x); names(x) = membernames; x })
}

satf_gridsearch <- function(start, constraints, ...) {
  rcpp_deinitialize_logLikFn()
  start.vec <- parameter_combinations(start)
  logLik.vec <- laply(start.vec, function(start) {
    constraints[names(start)] = start
    res = satf(start=start, constraints=c(constraints,start), ..., stepwise=F, .internal.init.optional=TRUE, .internal.cleanup=FALSE)
    reportifnot(!is.list(res), "All coefficient values have to be specified in a grid search.")
    res
  })
  rcpp_deinitialize_logLikFn()
  start.vec[[which.max(logLik.vec)]]
}

satf <- function(dv, signal, start, contrasts, bias, data, time, metric, trial.id=NULL, constraints=list(), 
                  optimize.incrementally=FALSE, reoptimize.dprime=TRUE, reoptimize.criterion=TRUE, reoptimize.corr=TRUE,
                  reoptimize.times=1, method="SUBPLEX", # "Nelder-Mead", #"BHHH",
                  debug=F, .likelihood.byrow=FALSE, .internal.init.optional=FALSE, .internal.cleanup=TRUE, doit=F)
{
  log = function(str) {
    if(debug) {
      cat(paste0(str,"\n"))
    }
  }
  log("---------- optimizing ---------")
  
  metric.permissible <- c('RMSD','R2','adjR2','logLik','logLikRaw')
  reportifnot(metric %in% metric.permissible, sprintf("'metric' has to be one of: %s", paste(metric.permissible, collapse=", ")))

  deinitialize_logLikFn <- function() T
  if(.internal.cleanup)
    deinitialize_logLikFn <- rcpp_deinitialize_logLikFn 
  
  # Check 'data' parameter
  reportifnot(nrow(data) > 0, "Parameter 'data' needs to consist of several rows.")
  
  satf.coefnames.core <- c('asymptote', 'invrate', 'intercept')
  bias.coefnames.core <- c('bias.max', 'bias.invrate', 'bias.intercept','bias.min')
  
  return_value = function(estimates, LL) {
    # if( is.na(LL) ) estimates = estimates*NaN
    list(estimates=estimates, LL=LL)#, SE=SE, ci.upper=ci.upper, ci.lower=ci.lower)
  }
  
  # set defaults for start values and constraints
  start = set_start_defaults(start, set.corr.mrsat=!is.null(trial.id))
  constraints = set_constraints_defaults(constraints)
  
  # check parameters
  params <- translate.parameters(data=data, dv=dv, contrasts=contrasts, bias=bias,
                                 signal=signal, time=time, trial.id=trial.id)

  skip.initialization <- .internal.init.optional && rcpp_is_initialized_logLikFn()
  
  # initialize start parameters and create constraint matrix if necessary
  if(!skip.initialization) {
    dm <- init_designmatrix(data=data, contrasts=params$contrasts, bias=params$bias, cnames=params$cnames,
                            satf.coefnames.core=satf.coefnames.core, bias.coefnames.core=bias.coefnames.core)
    coefnames <- colnames(dm$dm)
  } else {
    coefnames <- rcpp_get_coef_names()
  }
  
  # initialize start parameters and create constraint matrix  
  coefs <- init_coefs_and_constraints(coefnames=coefnames, start=start, constraints=constraints,
                                      coreparams=c(satf.coefnames.core, bias.coefnames.core))
#  # TODO: remove
#  coefs$constraints[,'upper'] = Inf
#  coefs$constraints[,'lower'] = -Inf
#  print(coefs$constraints)
  
  # initialize the C++ optimization routine
  if( !skip.initialization ) {
    rcpp_initialize_logLikFn(params$dv, dm$dm, dm$dm.coef.cnt, coefs$constraints, data, params$cnames)
  } else {
    rcpp_update_constraints_logLikFn(coefs$constraints)
  }
  
  start <- rcpp_unconstrain_coefs( coefs$start )
  names(start) <- names( coefs$start )
  
  if(all(names(start) %in% coefs$fixed.coefs)) {
    if(.likelihood.byrow) {
      logLik = compute_logLikFn(start, TRUE, TRUE)
      deinitialize_logLikFn()
      return( return_value( rcpp_constrain_coefs(start), logLik) )
      
    } else {
      logLik = compute_logLikFn(start, FALSE, TRUE)
      deinitialize_logLikFn()
      return( return_value( rcpp_constrain_coefs(start), logLik) )
    }
  }

  # make sure the start values yield a valid likelihood
  logLik = compute_logLikFn(start)
  if( is.na(logLik) ) {
    if(.likelihood.byrow)
     	res = compute_logLikFn(start, TRUE, FALSE)
    start = rcpp_constrain_coefs(start)
    deinitialize_logLikFn()
    warning("Got NA on first iteration. Adjust start parameters.")
    return( return_value( start, logLik) )
    
  } else if( is.infinite(logLik) ) {
    if(.likelihood.byrow)
      logLik = compute_logLikFn(start, TRUE, FALSE)
    start = rcpp_constrain_coefs(start)
    deinitialize_logLikFn()
    warning("Got Inf or -Inf on first iteration. Adjust start parameters.")
    return( return_value( start, logLik) )
  }

# TODO: Re-enable and do it again.
#  if(doit ==T) {
#  print( compareDerivatives(f=compute_logLikFn, grad=compute_logLikFn_gradient, t0=start) )
#  stop()
#  }
  
  optimize_subset <- function(variable, start, all=TRUE, selection=c(), method="Nelder-Mead") {
    .optimize_subset(variable, coefs$fixed.coefs, method, start, debug, data, selection)
  }

  log_step_n = function(n) {
    if(debug) {
      cat("-------------------\n")
      cat(sprintf("optimization step %d\n", n))
      cat("-------------------\n")
    }
  }
  

  original.start = start

  n.step = 0
  if(optimize.incrementally)
  {
    # ignore correlation parameter for now, if specified
    rcpp_set_coef_values( c(corr.mrsat=0) )
    
    # TODO: Implement the following.
    reportifnot(!.internal.init.optional, "The package does not support optimize.incrementally=T with .internal.init.optional=T. It's on my to-do list.")
    coeforder = append(dm$params.criterion, dm$params.dprime)
    cur.start = start
    # optimize in steps on subsets of the data
    for(n.step in 1:length(coeforder)) {
      log_step_n(n.step)
      selection.variables = c( unlist(coeforder[1:n.step]),  'corr.mrsat' )
      cur.start = optimize_subset(variable=coeforder[[n.step]], start=cur.start, selection=selection.variables) 
      if( any(is.na(cur.start)) ) {
        break;
      }
    }
    if( !any(is.nan(cur.start)) ) {
        start = cur.start
    }
    # optimize all of them
    log_step_n(n.step+1)
    start = optimize_subset(variable=c(unlist(dm$params.dprime), unlist(dm$params.criterion)), start=start, method=method) 
    
    # re-enable and optimize the correlation coefficient
    rcpp_reset_coef_ranges( 'corr.mrsat' )
    log_step_n(n.step+2)
    estimates = optimize_subset(variable='corr.mrsat', start=start, method=method) 
    if(any(is.na(estimates))) # Optimizing subsets of parameters has failed. Optimize over all of them together.
    {
      start = original.start
      reoptimize.dprime = TRUE
      reoptimize.criterion = TRUE
      reoptimize.corr = TRUE
    } else {
      start = estimates
    }
  }

    
  # reoptimize dprime or more parameters
  free.variables = c()
  if(reoptimize.dprime)    free.variables = c(free.variables, unlist(dm$params.dprime))
  if(reoptimize.criterion) free.variables = c(free.variables, unlist(dm$params.criterion))
  if(reoptimize.corr)      free.variables = c(free.variables, 'corr.mrsat') 
  for(i in 1:reoptimize.times) {  
    log_step_n(n.step+i+2)
    estimates = optimize_subset(variable=free.variables, start=start, method=method)
    start = estimates
  }

  if(.likelihood.byrow) {
    logLik = compute_logLikFn(coefs=estimates, by_row=TRUE)
    estimates = rcpp_constrain_coefs(estimates)
    deinitialize_logLikFn()
    return( return_value(estimates, logLik) )
  } 

  # recompute the likelihood for the entire dataset
  logLik <- compute_logLikFn( coefs= estimates, by_row=FALSE)
  constrained.estimates = rcpp_constrain_coefs(estimates)
  deinitialize_logLikFn()

  return( return_value(estimates=constrained.estimates, logLik) )
}

.optimize_subset <- function(variable, fixed, method, start, debug, data, selection, nosummary=TRUE)
{
  # NOTE: 'fixed' overrides 'variable'
  # extract start values if necessary
  original.start = start
  if(is.list(start))
    start = coef(start)
  
  # determine which variables should vary, and which are fixed
  variable.coefnames = setdiff(variable, fixed)
  fixed.coefnames = setdiff(names(start), variable.coefnames)
  fixed = names(start) %in% fixed.coefnames
  
  # if all are fixed, return the original coefs
  if( all(fixed) )
    return(original.start)
  
  n.free = sum(!fixed)
  if(n.free == 1 && method=="Nelder-Mead") {
    method = "CG" # Nelder-Mead can't handle one-parameter optimization
  }
    
  fnLogLik = compute_logLikFn
  
  # set gradient, etc.
  if(method == "Nelder-Mead")
    fnLogLikGradient = NULL
  else if(method == "BHHH")
    fnLogLikGradient = compute_logLikFn_gradient_BHHH
  else
    fnLogLikGradient = compute_logLikFn_gradient
  
  if(debug) print.level=0
  else      print.level=0
  
  # select the required data subset
  if(length(selection) > 0) {
    rcpp_select_coef_subset( selection )
  }

  if(debug) {
    variable.coefnames <- variable.coefnames[variable.coefnames%in%names(start)]
    cat(sprintf("optimizing: %s\n", paste(variable.coefnames, collapse=', ') ))
    cat(sprintf("data points: %d\n", length(rcpp_return_selection()) ))
  }
  
  generate_parscale <- function(start) {pmax(abs(start), .1)}
  
  logLik = compute_logLikFn(start)
  if(is.na(logLik) || is.infinite(logLik)) {
    estimates = start*NaN     
    rcpp_reset_selection()
    return(estimates)
  }
  if(method=="SUBPLEX") {
    reportifnot(all(!fixed), "Fixed parameters not yet supported when using SUBPLEX.")
    res = nloptr(x0=start, eval_f=fnLogLik, by_row=FALSE, force_update=FALSE, invert=TRUE, opts=list(algorithm="NLOPT_LN_SBPLX", maxeval=10^6))
    # TODO: translate status to code (maxLik codes) if necessary
    mapped.res = list(code=res$status, iterations = res$iterations, estimates=res$solution) 
    names(mapped.res$estimates) = names(start)
    res = mapped.res
    nosummary = TRUE
    coef = function(obj) obj$estimates
  } else if(method=="Nelder-Mead") {
    res = maxLik(logLik=fnLogLik, grad=fnLogLikGradient, start=start, fixed=fixed,
                 iterlim=10^6, method=method, print.level=print.level, parscale=generate_parscale(start))
  } else {
    res = maxLik(logLik=fnLogLik, grad=fnLogLikGradient, start=start, fixed=fixed,
                 iterlim=10^6, method=method, print.level=print.level)
  }
  i = 0
  switched_to_NM = FALSE
  while(res$code %in% c(3)) {  ## 100: Initial value out of range. 3: Boundary of parameter space.
    warning("Optimization routine got code 3. Restarting with Nelder-Mead.")
    cat("Optimization routine got code 3. Restarting with Nelder-Mead.\n")
    method = "Nelder-Mead"
    if(!nosummary)
      print(summary(res))
    switched_to_NM = TRUE
    i = i + 1
    if(i > 10) {
      break;
    }
    start = coef(res)
    if(method=="Nelder-Mead") {
      res = maxLik(logLik=fnLogLik, grad=fnLogLikGradient, start=start, fixed=fixed,
                   iterlim=10^6, method=method, print.level=print.level, parscale=generate_parscale(start))
    } else {
      res = maxLik(logLik=fnLogLik, grad=fnLogLikGradient, start=start, fixed=fixed,
                   iterlim=10^6, method=method, print.level=print.level)
    }
    if(debug && !nosummary) {
      print(summary(res))
    }
  }
#  if(switched_to_NM) {
#    start = coef(res)
#    res = maxLik(logLik=fnLogLik, grad=fnLogLikGradient, start=start, fixed=fixed,
#                 iterlim=10^6, method=method, print.level=print.level)
#  }
  if(debug) {
    cat(sprintf("method: %s\n", method))
    cat(sprintf("code: %d\n", res$code))
    cat(sprintf("iterations: %d\n", res$iterations))
    if(!nosummary)
      print(summary(res))
  }
  if(res$code %in% c(3,100)) { ## 100: Initial value out of range. 3: Boundary of parameter space.
    estimates = start*NaN     
    se = start*NaN
  } else {
    estimates = coef(res)
    if(nosummary)
      se = start*NaN
    else
      se = coef(summary(res))[,'Std. error']
  }
  
  if(debug) {
    cat("coefs\n")
    constrained.coefs = rcpp_constrain_coefs( estimates )
    idx = which(names(constrained.coefs)%in%variable.coefnames)
    names(constrained.coefs)[idx] = paste0('*',names(constrained.coefs)[idx],'*')
    newnames = names(constrained.coefs)
    constrained.coefs = sprintf("%.2f", constrained.coefs)
    names(constrained.coefs) = newnames
    print(constrained.coefs)
    old.LL = compute_logLikFn( coefs=start, by_row=FALSE)
    new.LL = compute_logLikFn( coefs=estimates, by_row=FALSE)
    cat(sprintf("LL improved by %.2f (old LL = %.2f, new LL = %.2f)\n", new.LL-old.LL, old.LL, new.LL))
  }
  rcpp_reset_selection()
  estimates
}
