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
  
  pc222 <- list(asym = c(1, 2), rate = c(1, 2), incp = c(1, 2))
  fit222 <- fit.SATcurve(data, par.cond = pc222)
  sum_curve <- summary.SATcurve(fit222)
  
  if(show_plot) plot(fit222, main = "222")
  
  data %>% ungroup() %T>% 
    {.$true_intercept1 = intercept[1]} %T>% {.$true_intercept2 = intercept[2]} %T>%
    {.$true_rate1 = rate[1]} %T>% {.$true_rate2 = rate[2]} %T>%
    {.$true_asymptote1 = asymptote[1]} %T>% {.$true_asymptote2 = asymptote[2]} %T>%
    {.$intercept1 = sum_curve$incp1} %T>% {.$intercept2 <- sum_curve$incp2} %T>%
    {.$rate1 = sum_curve$rate1} %T>% {.$rate2 = sum_curve$rate2} %T>%
    {.$asymptote1 = sum_curve$asym1} %T>% {.$asymptote2 = sum_curve$asym2} %T>%
    {.$R2 = sum_curve$R2} %T>% {.$adjR2 = sum_curve$adjR2} %T>%
    {.$AIC = -1*sum_curve$AIC} %>%
    dplyr::select(true_intercept1, true_intercept2, true_rate1, true_rate2, true_asymptote1, true_asymptote2,
                  R2, adjR2, AIC, intercept1, intercept2, rate1, rate2, asymptote1, asymptote2)
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

n_simulations <- 20
n_participants <- 15
n_trials_per_interval <- 50

time <- seq(0, 5.6, 0.35)

avg_incp <- 0.4
avg_rate <- 1
avg_asymp <- 3

delta_intercept = 0.05
delta_rate = 0.05
delta_asymptote = 0.05

intercept <- avg_incp + c(-.5, .5)*delta_intercept
rate <- avg_rate + c(-.5, .5)*delta_rate
asymptote <- avg_asymp + c(-.5, .5)*delta_asymptote

# For now, we assume that all participants are the same. We'll relax this assumption
# later. However, let's prepare the code for between-participants differences in 
# SATF parameters: Let's keep vectors of by-participant parameters.
intercepts <- rep(list(intercept), n_participants)
rates <- rep(list(rate), n_participants)
asymptotes <- rep(list(asymptote), n_participants)

# Run one simulated experiment with participant i
i <- 1
estimates <- sim_participant(n_trials_per_interval, time,
                             intercepts[[i]], rates[[i]], asymptotes[[i]], 
                             show_plot = FALSE)
estimates
# TODO: The above call should return *one* row. Fix. 


### FYI: (Below are the most elegant ways to loop when each iteration returns a
###       data frame, and if they are supposed to be bound by row.)
###
### How to loop over participants
tmp <- plyr::ldply(1:n_participants, function(i) {
    # ... sim_participant ...
    tmp <- data.frame(a = 1, b = 2, c = 3)
    tmp$participant_id <- i
})

### How to loop over simulations and participants
tmp <- plyr::ldply(1:n_simulations, function(j) {
  tmp2 <- plyr::ldply(1:n_participants, function(i) {
    # ... sim_participant ...
    tmp <- data.frame(a = 1, b = 2, c = 3)
    tmp$participant_id <- i
  })
  tmp2$simulation_id <- j
  tmp2
})

View(estimates)