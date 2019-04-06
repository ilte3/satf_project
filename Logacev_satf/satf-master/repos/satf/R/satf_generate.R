# TODO: remove the colums 'dprime', 'criterion', 'noise', and 'dprime.cur' from generated data frames,
#       or at least give them better names.

process.arg <- function(arg, time) {
  if(typeof(arg) == "closure") arg = arg(time)
  if(length(arg) == 1 ) arg = rep(arg, length(time))
  stopifnot(length(arg) == length(time) )
  arg
}


satf_generate_condition <- function(dprime, criterion, time, n.per.time, label, rho=0)
{
  if(rho==0)
    trial.ids = 1:(n.per.time*length(time))
  
  else
    trial.ids = rep(1:n.per.time, each=length(time))

  intervals <- 1:length(time)
  data <-   data.frame( condition=label, interval=intervals, time=time, 
                        trial.id=trial.ids, dprime=dprime, criterion=criterion )
  data$cur.epsilon <- rnorm( nrow(data) )
  if(rho != 0)
  data$cur.epsilon = rcpp_correlate(data$trial.id, data$cur.epsilon, rho)
  data$cur.psi <- data$dprime + data$cur.epsilon
  data$response <- data$cur.psi > data$criterion
  data[,c('condition', 'interval', 'time', 'trial.id', 'response')]
}

# 'criterion' is what Wickens calls 'lambda', and 'bias' is what he calls 'lambda_center'

satf_generate <- function(dprime, criterion=NULL, bias=NULL, time, n.per.time, rho=0, label="condition1") 
{
  if(is.null(criterion) && is.null(bias) )
    stop("bias or criterion have to be provided.")
  
  else if(is.null(criterion) && is.null(bias) )
    stop("Only one of bias and criterion may be provided.")
  
  dprime <- process.arg(dprime, time)
  if(!is.null(criterion)) {
    criterion <- process.arg(criterion, time)
  } else  {
    bias = process.arg(bias, time)
    criterion = bias + dprime/2
  }
  data0 <- satf_generate_condition(criterion=criterion, dprime=0, time=time, 
                                   n.per.time=n.per.time, label=label, rho=rho)
  data0$signal <- 0
  data1 <- satf_generate_condition(criterion=criterion, dprime=dprime,  time=time, 
                                   n.per.time=n.per.time, label=label, rho=rho)
  data1$signal <- 1
  data1$trial.id = data1$trial.id + max(data0$trial.id)
  rbind(data0, data1)
}
