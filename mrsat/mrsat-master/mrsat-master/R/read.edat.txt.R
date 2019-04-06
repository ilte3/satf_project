read.edat.txt <- function(file, collapse=FALSE)
{
  out <- list()
  for(fp in file)
  {
    line1 <- scan(fp, what="character", sep="\n",strip.white=TRUE, nlines = 1, quiet= TRUE)
    skip <- ifelse((grepl("emrg$", line1) || grepl("emrg$", line1)), 1, 0)
    out[[fp]] <- read.delim(fp,skip=skip)
    
  }
  if(collapse || length(file) == 1)
  {
    out <- do.call(rbind, out)
    rownames(out) <- NULL
  }
  out
}

read.param <- function(file, quiet = TRUE)
{
	if(!file.exists(file)) stop('parfile not found')
	npar <- scan(file = file, what = integer(0), nmax = 1, quiet=quiet)
	start <- scan(file = file, what = numeric(0), nmax = npar, sep = ",", skip = 1, quiet=quiet)
	lower <- scan(file = file, what = numeric(0), nmax = npar, sep = ",", skip = 2, quiet=quiet)
	upper <- scan(file = file, what = numeric(0), nmax = npar, sep = ",", skip = 3, quiet=quiet)
	param <- data.frame(start=start, lower=lower, upper=upper)
  param
}


read.dat <- function(file, quiet = TRUE)
{
    if(!file.exists(file)) stop('datfile not found')
    npoints <- scan(file = file, what = integer(0), nmax=1, quiet=quiet)
    
    d <- read.table(file,skip =1)
    dat <- matrix(c(t(d)), nrow=npoints, 
        dimnames = list(NULL, c("lags", "dprimes", "condition")))
    as.data.frame(dat)
}

get.param <- function(par.cond, asym=c(1, .2, 4.7), rate=c(3, .1, 12), incp=c(.25, .1, 3.0), random.start = FALSE, auto.asym = FALSE, data)
{
  if(is.list(par.cond)){
    if(any(!(c("asym", "rate", "incp") %in% names(par.cond)))) stop("some element missing in par.cond")
    rep.n <- sapply(par.cond, function(x){length(unique(x))})
  } 
  if(is.matrix(par.cond)) rep.n <- apply(par.cond, 2, function(x){length(unique(x))})

  param <- cbind(asym, rate, incp)
  if(any(param[1,] <= param[2,] || param[1,] >= param[3,])) stop("start should be within lower and upper")
  param <- data.frame(start=rep(param[1,], times=rep.n), lower=rep(param[2,], times=rep.n), upper=rep(param[3,], times=rep.n))
  if(random.start == TRUE) param$start <- runif(nrow(param), min=param$lower, max=param$upper)
  
  if(auto.asym){
    if(missing(data)) stop('needs data for auto.asym to work')
    if(!inherits(data, "mrsat.data")) stop("data needs to be class 'mrsat.dat' for auto.asym to work")
    n.asym <- length(unique(par.cond$asym))
    gt.bin <- max(data$bin) - 4 + 1
    cond <- as.numeric(factor(data$condition))
    grp <- par.cond$asym[cond]
    fix.asym <- tapply(data$dprimes[data$bin>gt.bin], grp[data$bin>gt.bin], mean)
    param$start[1:n.asym] <- fix.asym
  }
  param
}
