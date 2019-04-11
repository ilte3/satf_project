library(tidyverse)
library(magrittr)

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

satf_gen_cond <- function(mu, criterion, time, n, label) {
  stopifnot(length(mu) == length(criterion) && length(criterion) == length(time))
  trial = rep(1:n, each = length(time))
  interval <- 1:length(time)
  df <- data.frame(condition = label, interval = interval, time = time, trial = trial, mu = mu, criterion = criterion)
  df$p_yes <- p_yes(mu, criterion)
  df$response <- sapply(df$p_yes, function(p_yes) r_yes(n = 1, p_yes))
  df[, c('condition', 'interval', 'time', 'trial', 'response')]
}

satf_gen <- function(time, n, intercept, rate, asymptote, label = "condition1") {
  dprime <- satf(time, intercept, rate, asymptote)
  criterion <- 0.5*dprime
  df0 <- satf_gen_cond(mu = dprime*0, criterion = criterion, time = time, n = n, label = label)
  df1 <- satf_gen_cond(mu = dprime, criterion = criterion, time = time, n = n, label = label)
  df0$is_signal <- 0
  df1$is_signal <- 1
  df1$trial = df1$trial + max(df0$trial)
  rbind(df0, df1)
}