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

test_gen <- satf_gen(time = seq(0, 4.5, 0.5), n = 1000, intercept = 0.3, rate = 0.8, asymptote = 1, label = "condition1")
test_gen2 <- satf_gen(time = seq(0, 4.5, 0.5), n = 1000, intercept = 0.4, rate = 0.9, asymptote = 2, label = "condition2")
test_gen3 <- satf_gen(time = seq(0, 4.5, 0.5), n = 1000, intercept = 0.5, rate = 1, asymptote = 3, label = "condition3")
test_gen_bound <- rbind(test_gen, test_gen2, test_gen3)

satf_calc <- function(data) {
  calc_hit <- data.frame(data %>% filter(is_signal == 1)) %>% group_by(condition, interval) %>% 
    dplyr::summarize(hit = sum(response == TRUE))
  calc_miss <- data.frame(data %>% filter(is_signal == 1)) %>% group_by(condition, interval) %>% 
    dplyr::summarize(miss = sum(response == FALSE))
  calc_falarm <- data.frame(data %>% filter(is_signal == 0)) %>% group_by(condition, interval) %>% 
    dplyr::summarize(falarm = sum(response == TRUE))
  calc_creject <- data.frame(data %>% filter(is_signal == 0)) %>% group_by(condition, interval) %>% 
    dplyr::summarize(creject = sum(response == FALSE))
  calc <- left_join(calc_hit, calc_miss) %>% left_join(calc_falarm) %>% left_join(calc_creject)
  colnames(calc) <- c("condition", "interval", "hit", "miss", "falarm", "creject")
  calc
}

test_calc <- satf_calc(data = test_gen_bound)
View(test_calc)

cmp_dprime <- function(data, hit, falarm, miss, creject) {
  n_signal <- hit + miss
  n_noise <- falarm + creject
  p_hit <- hit/n_signal
  p_falarm <- falarm/n_noise
  data$dprime <- qnorm(p_hit) - qnorm(p_falarm)
  data$criterion <- -0.5*(qnorm(p_hit) + qnorm(p_falarm))
  data
}

test_dprimes <- cmp_dprime(data = test_calc, hit = test_calc$hit, miss = test_calc$miss,
                           falarm = test_calc$falarm, creject = test_calc$creject)
View(test_dprimes)