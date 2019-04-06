satf <- function(t, intercept, rate, asymptote) {
  ifelse(t > intercept, asymptote*(1-exp(-rate*(t-intercept))), 0)
}

p_yes <- function(mu, criterion) {
  1 - pnorm(criterion, mean = mu, sd = 1)
}

r_yes <- function(n, p_yes) {
  stopifnot(length(n) == 1 && length(p_yes) == 1)
  runif(n) <= p_yes
}

### test arguments
mu = 3
criterion = 0.5*mu
time = 3.5
n = 10

(with(df[1:nrow(df),], {
p_yes (mu, criterion)
}))

x <- p_yes (3, 1.5)

r_yes(1, x)

?pnorm

pnorm(180, 172, 10)
pnorm(10, 130, 30)

pnorm(c(180, 10), c(172, 130), c(10, 30))

sapply(df$p_yes, function(p_yes) r_yes(n = 1, p_yes))

###

satf_gen_cond <- function(mu, criterion, time, n) {
  stopifnot(length(mu) == length(criterion) && length(criterion) == length(time))
  trial = rep(1:n, each = length(time))
  interval <- 1:length(time)
  df <- data.frame(interval = interval, time = time, trial = trial, mu = mu, criterion = criterion)
  df$p_yes <- p_yes(mu, criterion)
  df$response <- sapply(df$p_yes, function(p_yes) r_yes(n = 1, p_yes))
  df[, c('interval', 'time', 'trial', 'response')]
  }

satf_gen_cond(seq(0, 9, 1), seq(0, 4.5, 0.5), seq(0, 9, 1), n = 10^4)

satf_gen <- function(time, n, intercept, rate, asymptote) {
  dprime <- satf(time, intercept, rate, asymptote)
  criterion <- 0.5*dprime
  df0 <- satf_gen_cond(mu = dprime*0, criterion = criterion, time = time, n = n)
  df1 <- satf_gen_cond(mu = dprime, criterion = criterion, time = time, n = n)
  df0$is_signal <- 0
  df1$is_signal <- 1
  df1$trial = df1$trial + max(df0$trial)
  rbind(df0, df1)
  }

satf_gen 

dprime <- function(t) satf(seq(0, 4.5, 0.5), 0.5, 1, 5)
n <- 100
time = seq(0, 4.5, 0.5)
intercept <- 0.5
rate = 1
asymptote = 3

my_df <- satf_gen(time = time, n = n, intercept = intercept, rate = rate, asymptote = asymptote)
View(my_df)

# ###
# iterate over intervals
# compute the number of hits, falarms, misses, correjects
# group_by summarize mutate - - - tidyr spread
# avoid high asymptotes
# ###

cmp_dprime <- function(hit, falarm, miss, creject) {
  n_signal <- hit + miss
  n_noise <- falarm + creject
  p_hit <- hit/n_signal
  p_falarm <- falarm/n_noise
  dprime <- qnorm(p_hit) - qnorm(1 - p_falarm)
  criterion <- 0.5*(qnorm(p_hit) + qnorm(1 - p_falarm))
  c(dprime = dprime, criterion = criterion)
}
