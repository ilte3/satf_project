
add.params <- function(formula, letter, startpars, start=0) {
  n <- n.terms(formula)
  if(n > 0) {
    names <- paste(letter, formula.terms(formula), sep='.')
    vec <- ifelse(!is.na(startpars[names]), startpars[names], start)
    names(vec) <- names
    return(vec)
  }
  return(c())
}

init_designmatrix <- function(data, contrasts, bias, cnames, satf.coefnames.core, bias.coefnames.core) {

  reportifnot(all(satf.coefnames.core %in% names(contrasts)), 
                  "Parameter 'contrasts' needs to contain 'asymptote', 'invrate', and 'intercept'.")
  if( !is.list(bias) ) {
    bias.contrast <- bias
    bias <- list()
    for(coreparam.bias in bias.coefnames.core)
      bias[[coreparam.bias]] = bias.contrast

  } else if( !all(bias.coefnames.core %in% names(bias)) ) {
    stop("Parameter 'bias' needs to contain formulae for 'bias.min', 'bias.max', 'bias.invrate', and 'bias.intercept', or to be a formula.")
    
  }

  # rearange start parameters, make sure the coefficients for every parameter are contiguous, 
  # and create the design matrix

  process_contrasts <- function(contrasts, default.colname) {
    dm.colnames = c();
    allparamnames = c();
    paramnames = list();
    for(paramname in names(contrasts)) {
      allparamnames = c(allparamnames, paramname)
      paramnames[['core']] = c(paramnames[['core']], paramname)
      dm.colnames = c(dm.colnames, default.colname)
      for(cur.colname in contrasts[[paramname]]) {
        cur.paramname = paste(paramname, cur.colname, sep='.')
        allparamnames = c(allparamnames, cur.paramname)
        paramnames[[cur.colname]] = c(paramnames[[cur.colname]], cur.paramname)
        dm.colnames = c(dm.colnames, cur.colname)
      }
    }
    list(paramnames=paramnames, allparamnames=allparamnames, dm.colnames=dm.colnames)
  }
  
  params.dprime = process_contrasts(contrasts, cnames[['signal']])
  params.criterion = process_contrasts(bias, '1')
  dm.coef.cnt = c(unlist(lapply(contrasts, length))+1, unlist(lapply(bias, length))+1)

  coef.names = c(params.dprime$allparamnames, params.criterion$allparamnames)
  dm.colnames = c(params.dprime$dm.colnames, params.criterion$dm.colnames)
  
  # create a column consisting of only 1's, for use with the bias
  data$'1' <- 1
  
  missing.columns = dm.colnames[!(dm.colnames %in% colnames(data))]
  if(length(missing.columns) > 0) {
    stop(sprintf("Missing columns <%s>", paste0(missing.columns, sep=",", collapse=",")))
  }
  
  # create the design matrix
  dm <- data[,dm.colnames]
  colnames(dm) <- coef.names
  
  # set all satf coefficients to 0 for noise trials
  satf.ncol <- sum(dm.coef.cnt[satf.coefnames.core])
  dm[,1:satf.ncol] <- dm[,1:satf.ncol]*data[, cnames['signal'] ]
  
  # add a corr.mrsat parameter if necessary
  if( 'trial.id' %in% names(cnames) ) {
    coef.names <- c(coef.names, "corr.mrsat")
    dm.coef.cnt <- c(dm.coef.cnt, corr.mrsat=1)
    trial.id <- data[[ cnames[['trial.id']] ]]
    dm$corr.mrsat <-  as.integer(c(F, trial.id[-1] == trial.id[-length(trial.id)]))
  }
  
  list(dm=as.matrix(dm), dm.coef.cnt=dm.coef.cnt, 
       params.dprime=params.dprime$paramnames, params.criterion=params.criterion$paramnames)
}



init_coefs_and_constraints <- function(coefnames, start, constraints, coreparams)
{
  reportifnot(is.vector(start) && !is.list(start), "Parameter 'start' needs to be a vector.")
  reportifnot(!is.null(start) || !any(is.na(start)), "Parameter 'start' was not provided or containts NAs.")

  # make sure start values were provided for all the necessary parameters
  missing.coreparams = coreparams[!(coreparams %in% names(start))]
  missing.coreparams.str = paste(missing.coreparams, collapse=", ")
  reportifnot( length(missing.coreparams) == 0, sprintf("Parameter 'start' does not contain start values for: %s", 
                                                        missing.coreparams.str))

  # warn about unused start values
  unused.startvalues = names(start[!(names(start) %in% coefnames)])
  unused.startvalues.str = paste(unused.startvalues, collapse=", ")
  if( length(unused.startvalues) > 0) {
      warning(sprintf("Parameter 'start' contains unused start values for: %s", unused.startvalues.str))
  }

  start = start[coefnames]
  names(start) = coefnames

  # construct the constraint matrix with a reversed the constraint list to make user-specified constraints take precedence over defaults
  # (the reason there can be several specifications for one parameter is that constraints names are specified by regex expressions)
  constraint.matrix = init_constraint_matrix(coefnames, rev(constraints))
  
  new_start_value <- function(constraint, oldvalue=NA) 
  {
    correction.fraction = .1
    # find the value closest to the old value within the boundaries (or to 0 is the old value was NA)
    if( is.na(oldvalue)) {
	if( constraint[1] <= 0 && constraint[2] >= 0) {
		return(0);
	} else {
		return( mean(constraint) );
}
    }

    if( all(is.infinite(constraint)) ) {
      newvalue = 0

    } else if( is.infinite(constraint[1]) ) {
      newvalue = constraint[2]-correction.fraction

    }  else if( is.infinite(constraint[2]) ) {
      newvalue = constraint[1]+correction.fraction

    } else if( constraint[1]==constraint[2] ) {
      newvalue = constraint[1]

    } else {
      if( abs(constraint[1]-oldvalue) > abs(constraint[2]-oldvalue) ) {
        newvalue = constraint[2]-correction.fraction
      } else {
        newvalue = constraint[1]+correction.fraction
      }
    }
    newvalue
  }

  # fill the constraints matrix
  start.changed = rep(FALSE, length(start))
  names(start.changed) = names(start)
  for(coefname in coefnames)
  {
    cur.start = start[[coefname]]
    cur.constraint = constraint.matrix[coefname,]
    if( is.na(cur.start) ) {
      start[[coefname]] = new_start_value(cur.constraint)

    } else if( cur.start < cur.constraint[['lower']] || cur.start > cur.constraint[['upper']] ) 
    {
      start[[coefname]] = new_start_value(cur.constraint, cur.start)
      start.changed[[coefname]] = TRUE
    } 
  }

  if( any( start.changed )) {
    start.changed = names(start.changed[start.changed])
    start.changed = paste0(start.changed, '=', start[start.changed],  collapse=", ")
    warning(sprintf("Start values changed to <%s>.", start.changed))
  }

  fixed.coefs <- coefnames[constraint.matrix[,'upper'] == constraint.matrix[,'lower']]
  
  return(list(constraints=constraint.matrix, start=start, fixed.coefs=fixed.coefs));
}


  
init_constraint_matrix <- function(coefnames, constraints) 
{
    # define functions for finding coefficient constraints in the constraint list
    is_match <- function(expr, name) regexpr(paste0('^',expr,'$'), text=name) > 0
    find_coefnames <- function(expr, coefnames) {
      res = c()
      for(coefname in coefnames) {
        if(is_match(expr, coefname)) {
          res = c(res, coefname)
        }
      }
      return(res)
    }

    add_constraint <- function(cmatrix, coefname, constraint)
    {
      # use the most specific constraints specified as long as they are coherent,
      # if they are incoherent, use the ones first in the reversed list

      # set to fixed value if constraint is fixed
      if(length(constraint) == 1) {
        cmatrix[coefname,] = constraint

      # or set boundaries, unless they are already set to a fixed value
      } else if(cmatrix[coefname, 'lower'] != cmatrix[coefname, 'upper']) 
      {
        # increase lower boundary
        if( constraint[1] > cmatrix[coefname, 'lower'] )
          cmatrix[coefname, 'lower'] = constraint[1]

        # decrease upper boundary
        if( constraint[2] < cmatrix[coefname, 'upper'] )
          cmatrix[coefname, 'upper'] = constraint[2]

        reportifnot(cmatrix[coefname, 'lower'] <= cmatrix[coefname, 'upper'],
          sprintf("Upper boundary for parameter %s is below lower boundary (%f; %f).", 
                  coefname, cmatrix[coefname, 'lower'], cmatrix[coefname, 'upper']))
      }
      cmatrix
    }

    # initialize constraints matrix
    constraint.matrix = matrix(c(-Inf, Inf), nrow=length(coefnames), ncol=2, byrow=TRUE, 
                               dimnames=list(coefnames, c('lower', 'upper')))

    for(expr in names(constraints))
    {
      constraint = constraints[[expr]]
      for(coefname in find_coefnames(expr, coefnames)) {
        constraint.matrix = add_constraint(constraint.matrix, coefname, constraint)
      }
    }
    constraint.matrix
}


# translate parameters from formula notation to string notation,
# and check that all columns actually exist
# TODO: This function must be rather slow, optimize.

translate.parameters  <- function(data, dv, contrasts, bias, signal, time, trial.id) {
  
  cnames <- list()

  # translate parameter: contrasts
  for(name in names(contrasts)) 
    contrasts[[name]] = check.formula.for.colname(data, paste('contrasts', name, sep='$'), contrasts[[name]], n.cols=NA)

  # translate parameter: contrasts
  if(is.list(bias)) {
    for(name in names(bias)) 
      bias[[name]] = check.formula.for.colname(data, paste('bias', name, sep='$'), bias[[name]], n.cols=NA)

  } else {
      bias = check.formula.for.colname(data, 'bias', bias, n.cols=NA)

  }
  
  # TODO: Make sure that trials are not discontinuous.
  # check column names: signal, time, trial.id
  cnames$signal <- check.formula.for.colname(data, 'signal', signal)
  stopifnot.binary( data[,cnames$signal] )
  
  cnames$time <- check.formula.for.colname(data, 'time', time)
  
  if(!is.null(trial.id))
    cnames$trial.id <- check.formula.for.colname(data, 'trial.id', trial.id)
  
  # translate parameter: dv
  if('response' %in% names(dv)) {
    dv[['response']] <- check.formula.for.colname(data, 'dv$response', dv$response)
    stopifnot.binary( data[, dv[['response']] ] )
    
  } else if( all(c('n.responses.yes','n.responses') %in% names(dv)) ) {
    for(name in names(dv))
      dv[[name]] = check.formula.for.colname(data, paste('dv', name, sep='$'), dv[[name]])
    
  } else {
    stop("No such dv implemented.")
    
  }
  dv <- unlist(dv)

  # return
  list(dv=dv, contrasts=contrasts, bias=bias, cnames=unlist(cnames) )
}

set_start_defaults <- function(start, set.corr.mrsat=FALSE) {
  start = default(start, 'bias.min', -1)
  start = default(start, 'bias.max', 1)
  start = default(start, 'bias.invrate', 1)
  start = default(start, 'bias.intercept', 0)
  if(set.corr.mrsat)
    start = default(start, 'corr.mrsat', .85)
  start
}

set_constraints_defaults <- function(constraints) {
  constraints = default(constraints, 'asymptote', c(0,5))
  #constraints = default(constraints, 'asymptote.*', c(-5,5))
  constraints = default(constraints, 'invrate', c(0,Inf))
  constraints = default(constraints, 'intercept', c(0,Inf))
  constraints = default(constraints, 'bias.invrate', c(0,Inf))
  # Albers and Kallenberg recommend [1/sqrt(2), 1] for their approximation
  constraints = default(constraints, 'corr.mrsat', c(1/sqrt(2), 1))
  constraints
}
