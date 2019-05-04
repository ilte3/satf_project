source("../default.R")
library(mrsat)

satf_calc <- function(data)
{
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

cmp_dprime <- function(data, hit, fa, miss, creject) {
  n_signal <- hit + miss
  n_noise <- fa + creject
  p_hit <- hit/n_signal
  p_fa <- fa/n_noise
  data$dprime <- qnorm(p_hit) - qnorm(p_fa)
  data$criterion <- -0.5*(qnorm(p_hit) + qnorm(p_fa))
  data %<>% dplyr::select(interval, time, hit, miss, n_signal, fa, creject, n_noise, dprime, criterion, condition) 
  data
}


time <- seq(0, 4.5, 0.5)
n = 100

test_gen <- satf_gen(time = time, n = n, intercept = 0.3, rate = 0.8, asymptote = 1, label = "noint")
test_gen2 <- satf_gen(time = time, n = n, intercept = 0.4, rate = 0.9, asymptote = 2, label = "obrel")
test_gen3 <- satf_gen(time = time, n = n, intercept = 0.5, rate = 1, asymptote = 3, label = "obrelsub")
test_gen_bound <- rbind(test_gen, test_gen2, test_gen3)


test_calc <- satf_calc(data = test_gen_bound)
head(test_calc)

test_dprimes <- cmp_dprime(data = test_calc, hit = test_calc$hit, miss = test_calc$miss,
                           fa = test_calc$fa, creject = test_calc$creject)
head(test_dprimes)


mrsat_dprimes <- test_dprimes %T>%
                      {.$hit.correction = 'none'} %T>% 
                      {.$fa.correction = 'none'} %>% 
                      dplyr::select(bin=interval, hit, hit.correction,
                                    hit.denom=n_signal, fa, fa.correction, fa.denom=n_noise, 
                                    lags=time, dprimes=dprime, condition)

pc333 <- list(asym = c(1, 2, 3), rate = c(1, 2, 3), incp = c(1, 2, 3))
fit333 <- fit.SATcurve(mrsat_dprimes, par.cond = pc333)
plot(fit333, main = "333")

pc311 <- list(asym = c(1, 2, 3), rate = c(1, 1, 1), incp = c(1, 1, 1))
fit311 <- fit.SATcurve(mrsat_dprimes, par.cond = pc311)
plot(fit311, main="311")

SATsummary.list(list(fit333, fit311))

# create a table with experiment information (participants and such)
# draw one incp rate asymp - "simulated participant"
# one sample - two conditions - fitsatcurve 333 ->>>>> extract params
fit333$fit$par
# create a function for all this
# asymp mean 3 std 1
# rate mean 1 std 0.1 sqr
# incp mean 0.4 std 0.1 sqr

# needs to return INFO
# input avgs, deltas
# how much we are off by