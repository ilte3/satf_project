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
  #df$p_yes <- p_yes(mu, criterion)
  #df$response <- sapply(df$p_yes, function(p_yes) r_yes(n = 1, p_yes))
  df$evidence_position <- rnorm(nrow(df), mean = mu, sd = 1)
  df$response <- df$evidence_position > criterion
  df %>% dplyr::select(condition, interval, time, trial, response)
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

cmp_acc_stat <- function(data)
{
  calc <- data %>% group_by(condition, interval, is_signal) %>% 
                   dplyr::summarize(p_yes = mean(response), time = mean(time), n = length(response))
  
  calc_signal <- calc %>% filter(is_signal == 1) %>% dplyr::select(condition, interval, time, hit = p_yes, n_signal = n) %>% dplyr::mutate(miss = 1-hit)
  calc_noise <- calc %>% filter(is_signal == 0) %>% dplyr::select(condition, interval, time, fa = p_yes, n_noise = n) %>% dplyr::mutate(creject = 1-fa)
  
  left_join(calc_signal, calc_noise, by = c("condition", "interval", "time")) %>% 
            dplyr::select(interval, time, hit, miss, n_signal, fa, creject, n_noise, condition )
}

cmp_dprime <- function(data) {
  n_signal <- data$hit + data$miss
  n_noise <- data$fa + data$creject
  p_hit <- data$hit/n_signal
  p_fa <- data$fa/n_noise
  data$dprime <- qnorm(p_hit) - qnorm(p_fa)
  data$criterion <- -0.5*(qnorm(p_hit) + qnorm(p_fa))
  data %<>% dplyr::select(interval, time, hit, miss, n_signal, fa, creject, n_noise, dprime, criterion, condition) 
  data
}

mrsat_fitcurve <- function(data, pc = list(asym = c(1, 2), rate = c(1, 2), incp = c(1, 2)), show_plot = FALSE) {
  data <- data %T>% 
  {.$hit.correction = 'none'} %T>% 
  {.$fa.correction = 'none'} %>% 
    dplyr::select(bin = interval, hit, hit.correction,
                  hit.denom = n_signal, fa, fa.correction, fa.denom = n_noise, 
                  lags = time, dprimes = dprime, condition)
  
  fit <- fit.SATcurve(data, par.cond = pc)
  sum_curve <- summary.SATcurve(fit)
  
  if(show_plot) plot(fit, main = "222")
  
  n_unique <- function(x) length(unique(x))
  
  model_id <- sapply(pc, n_unique) %>% paste(collapse = "-")

  model_fit <- sum_curve[c("incp1", "incp2", "rate1", "rate2", "asym1", "asym2", "R2", "adjR2", "AIC")]
  model_fit$model <- model_id
  model_fit %<>% dplyr::rename(intercept1 = incp1, intercept2 = incp2, asymptote1 = asym1, asymptote2 = asym2)
  
  model_fit
}

sim_participant <- function(n, time, intercept, rate, asymptote, show_plot = FALSE) 
{
  responses1 <- satf_gen(time = time, n = n, intercept = intercept[1], rate = rate[1], asymptote = asymptote[1], label = "condition1")
  responses2 <- satf_gen(time = time, n = n, intercept = intercept[2], rate = rate[2], asymptote = asymptote[2], label = "condition2")
  responses <- rbind(responses1, responses2)
  
  acc_stat <- cmp_acc_stat(data = responses)
  acc_stat <- cmp_dprime(data = acc_stat)
  
  res111 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 1), rate = c(1, 1), incp = c(1, 1)), show_plot = FALSE)
  res112 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 1), rate = c(1, 1), incp = c(1, 2)), show_plot = FALSE)
  res121 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 1), rate = c(1, 2), incp = c(1, 1)), show_plot = FALSE)
  res122 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 1), rate = c(1, 2), incp = c(1, 2)), show_plot = FALSE)
  res211 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 2), rate = c(1, 1), incp = c(1, 1)), show_plot = FALSE)
  res221 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 2), rate = c(1, 2), incp = c(1, 1)), show_plot = FALSE)
  res212 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 2), rate = c(1, 1), incp = c(1, 2)), show_plot = FALSE)
  res222 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 2), rate = c(1, 2), incp = c(1, 2)), show_plot = FALSE)
  res <- dplyr::bind_rows(res111, res112, res121, res122, res211, res221, res212, res222)
  
  res %<>% cbind(data.frame(true_intercept1 = intercept[1], true_intercept2 = intercept[2],
                           true_rate1 = rate[1], true_rate2 = rate[2],
                           true_asymptote1 = asymptote[1], true_asymptote2 = asymptote[2]))
  
  res %>% dplyr::select(model, true_intercept1, true_intercept2, true_rate1, true_rate2,
                        true_asymptote1, true_asymptote2, R2, adjR2, AIC,
                        intercept1, intercept2, 
                        rate1, rate2, 
                        asymptote1, asymptote2)
}
