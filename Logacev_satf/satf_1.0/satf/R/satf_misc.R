
# TODO: Use as.quoted() for formulas.
check.formula.for.colname <- function(data, name, formula, n.cols=1) {
  reportifnot( is.formula(formula), sprintf("Parameter '%s' needs to be a formula.", name) )
  colname <- formula.terms(formula)
  if(!is.na(n.cols))
    reportifnot( length(colname) == n.cols, sprintf("Parameter '%s' needs specify exactly %d column.", name, n.cols) )
  if(length(colname))
    check.columns( colname, data )
  colname
}

get.column <- function(data, name, formula, as.type) {
  colname = check.formula.for.colname(data, name, formula)
  as.type( data[, colname ] )
}

stopifnot.binary <- function(x) stopifnot( all(x %in% c(0,1)) )
stopifnot.binary.na <- function(x) stopifnot( all(x %in% c(0,1,NA)) )


reportifnot <- function(cond, text) {
  if(!cond)
    stop(text)
}

check.columns <- function(columns, data) {
  present <-  columns %in% colnames(data)
  if(!all(present)) {
    missing <- columns[!present]
    missing <- paste("'",missing,"'",sep="", collapse=", ")
    msg <- paste("Missing columns in data:", missing)
    stop(msg)
  }
}

default <- function(lst, name, val) {
  if(!name %in% names(lst)) lst[[name]] <- val; 
  lst
}

optim.to.precision <- function(start, optim.digits, control, method="Nelder-Mead", ...)
{
  run.optim <- function(start) optim(par=start, method=method, control=control, ...)
  res <- run.optim(start)
  if(method == "SANN")
    return(res)
  old.cnt <- NULL
  while(TRUE) {
    old.value <- res$value
    old.cnt <- paste(old.cnt, res$counts[['function']], sep=' ')
    res <- run.optim(res$par)
    if(is.na(optim.digits))
      return(res)
    if(round(old.value, optim.digits) == round(res$value, optim.digits)) {
      res$counts[['function']]  <- paste(old.cnt, res$counts[['function']], sep=' ')
      return(res)
    }
  }
}

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

formula.terms <- function(formula) { 
  terms.formula <- terms(formula)
  stopifnot(attr(terms.formula, "intercept")==1)
  attr(terms.formula, "term.labels")
}

n.terms <- function(formula) length(formula.terms(formula))

SATF <- function(t, asymptote, invrate, intercept)
  (t >= intercept)*(asymptote*(1-exp(-1/invrate*(t-intercept)))) 

