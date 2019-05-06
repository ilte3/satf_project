source("../default.R")
library(mrsat)

cmp_acc_stat <- function(data) {
  calc_hit <- data %>% filter(is_signal == 1) %>% group_by(condition, interval) %>% 
    dplyr::summarize(hit = mean(response), time = mean(time), n_signal = length(response))
  
  calc_miss <- data %>% filter(is_signal == 1) %>% group_by(condition, interval) %>% 
    dplyr::summarize(miss = mean(!response), time = mean(time))
  
  calc_fa <- data %>% filter(is_signal == 0) %>% group_by(condition, interval) %>% 
    dplyr::summarize(fa = mean(response), time = mean(time), n_noise = length(response))
  
  calc_creject <- data %>% filter(is_signal == 0) %>% group_by(condition, interval) %>% 
    dplyr::summarize(creject = mean(!response), time = mean(time))
  
  calc <- left_join(calc_hit, calc_miss) %>% left_join(calc_fa) %>% left_join(calc_creject)
  calc %<>% dplyr::select(interval, time, hit, miss, n_signal, fa, creject, n_noise, condition )
  calc
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

mrsat_fitcurve <- function(data, show_plot = FALSE) {
  data <- data %T>% 
    {.$hit.correction = 'none'} %T>% 
    {.$fa.correction = 'none'} %>% 
    dplyr::select(bin = interval, hit, hit.correction,
                  hit.denom = n_signal, fa, fa.correction, fa.denom = n_noise, 
                  lags = time, dprimes = dprime, condition)
  
  pc112 <- list(asym = c(1, 1), rate = c(1, 1), incp = c(1, 2))
  fit112 <- fit.SATcurve(data, par.cond = pc112)
  
  if(show_plot) plot(fit112, main = "112")
  
  asymptotes <- fit112$fit$par[1]
  rates <- fit112$fit$par[2]
  intercepts <- fit112$fit$par[3:4]
  
  data.frame(asymptotes, rates, intercepts)
}

sim_participant <- function(n, time, intercept, rate, asymptote, show_plot = FALSE) 
{
  responses1 <- satf_gen(time = time, n = n, intercept = intercept[1], rate = rate[1], asymptote = asymptote[1], label = "condition1")
  responses2 <- satf_gen(time = time, n = n, intercept = intercept[2], rate = rate[2], asymptote = asymptote[2], label = "condition2")
  responses <- rbind(responses1, responses2)
  
  acc_stat <- cmp_acc_stat(data = responses)
  acc_stat <- cmp_dprime(data = acc_stat)
  
  mrsat_fitcurve(acc_stat, show_plot = show_plot)
}

###################################################################################################################################

n <- 50
time <- seq(0, 5.6, 0.35)

avg_incp <- rnorm(1, 0.4, 0.1)^2
avg_rate <- rnorm(1, 1, 0.1)^2
avg_asymp <- rnorm(1, 3, 1)

delta_intercept = 0
delta_rate = 0
delta_asymptote = 0

intercept <- avg_incp + c(-.5, .5)*delta_intercept
rate <- avg_rate + c(-.5, .5)*delta_rate
asymptote <- avg_asymp + c(-.5, .5)*delta_asymptote

estimates <- sim_participant(n, time, intercept, rate, asymptote, show_plot = TRUE)
estimates

###

delta_intercept1 = 0
intercept1 <- avg_incp + c(-.5, .5)*delta_intercept1
estimates1 <- sim_participant(n, time, intercept1, rate, asymptote, show_plot = TRUE)

delta_intercept2 = 0.05
intercept2 <- avg_incp + c(-.5, .5)*delta_intercept2
estimates2 <- sim_participant(n, time, intercept2, rate, asymptote, show_plot = TRUE)

delta_intercept3 = 0.1
intercept3 <- avg_incp + c(-.5, .5)*delta_intercept3
estimates3 <- sim_participant(n, time, intercept3, rate, asymptote, show_plot = TRUE)