satf <- function(t, intercept, rate, asymptote) {
  ifelse(t > intercept, asymptote*(1-exp(-rate*(t-intercept))), 0)
}

satf_gen_cond <- function(dprime, criterion, time, n) {
  trial = rep(1:n, each = length(time))
  interval <- 1:length(time)
  df <- data.frame(interval = interval, time = time, trial = trial, dprime = dprime, criterion = criterion)
  p_yes_dprime <- function(dprime, criterion, is_signal) {
    1 - pnorm(crit, mean = is_signal*dprime, sd = 1)
  }
  r_yes <- function(n, p_yes) {
    stopifnot(length(n) == 1 && length(p_yes) == 1)
    runif(n) <= p_yes
  }
  df$response <- df$dprime > df$criterion
  df[, c('interval', 'time', 'trial', 'response')]
}

satf_gen <- function(time, n, intercept, rate, asymptote) {
  dprime <- satf(time, intercept, rate, asymptote)
  criterion <- 0.5*dprime
  df0 <- satf_gen_cond(dprime = 0, criterion = criterion, time = time, n = n)
  df1 <- satf_gen_cond(dprime = dprime, criterion = criterion, time = time, n = n)
  df0$is_signal <- 0
  df1$is_signal <- 1
  df1$trial = df1$trial + max(df0$trial)
  rbind(df0, df1)
}


dprime <- function(t) satf(0, 0.5, 1, 5)
n <- 10
time = seq(0, 4.5, 0.5)

my_df <- satf_gen(criterion = criterion, dprime = dprime, time=time, n=n)
View(my_df)

p_yes_dprime <- function(dprime, crit, is_signal) {
  1 - pnorm(crit, mean = is_signal*dprime, sd = 1)
}

dprime <- satf(seq(0, 4.5, 0.5), 0.5, 1, 5)
dprime

p_yes_dprime(dprime, 0.5*dprime, 0)
p_yes_dprime(dprime, 0.5*dprime, 1)

r_yes <- function(n, p_yes) {
  stopifnot(length(n) == 1 && length(p_yes) == 1)
  runif(n) <= p_yes
} 

x <- r_yes(20000, .3)
mean(x)

# simple alternative to a loop
sapply(1:10, function(x) x)
sapply(1:10, function(x) x/2)

sapply(seq(0,1,.1), function(p) r_yes(1, p))

sapply(seq(0, 4.5, 0.5), function(dprime) r_yes(10, dprime))

sapply(dprime, r_yes())

cmp_dprime <- function(hit, falarm, miss, creject) {
  n_signal <- hit + miss
  n_noise <- falarm + creject
  p_hit <- hit/n_signal
  p_falarm <- falarm/n_noise
  dprime <- qnorm(p_hit) - qnorm(1 - falarm)
  criterion <- 0.5*(qnorm(p_hit) + qnorm(1 - p_falarm))
  c(dprime = dprime, criterion = criterion)
}

r_yes_dprime <- for(value in dprime[1:10]) {
  res_dprime <- r_yes(10, value)
  print(res_dprime)
}


p_yes_dprime <- function(dprime, crit, is_signal) {
  1 - pnorm(crit, mean = is_signal*dprime, sd = 1)
}

r_yes <- function(n, p_yes) {
  stopifnot(length(n) == 1 && length(p_yes) == 1)
  runif(n) <= p_yes
} 

for(value in dprime[1:10]) {
  res_dprime <- r_yes(10, value)
}
