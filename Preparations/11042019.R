source("../default.R")

time <- seq(0, 4.5, 0.5)
n = 2000

test_gen <- satf_gen(time = time, n = n, intercept = 0.3, rate = 0.8, asymptote = 1, label = "noint")
test_gen2 <- satf_gen(time = time, n = n, intercept = 0.4, rate = 0.9, asymptote = 2, label = "obrel")
test_gen3 <- satf_gen(time = time, n = n, intercept = 0.5, rate = 1, asymptote = 3, label = "obrelsub")
test_gen_bound <- rbind(test_gen, test_gen2, test_gen3)

satf_calc <- function(data) {
  data <- rename(data, bin = interval, lags = time)
  calc_hit <- data %>% filter(is_signal == 1) %>% group_by(condition, bin) %>% 
              dplyr::summarize(hit = sum(response == TRUE)/sum(is_signal == 1), lags = mean(lags))
  calc_miss <- data %>% filter(is_signal == 1) %>% group_by(condition, bin) %>% 
               dplyr::summarize(miss = sum(response == FALSE)/sum(is_signal == 1), lags = mean(lags))
  calc_fa <- data %>% filter(is_signal == 0) %>% group_by(condition, bin) %>% 
             dplyr::summarize(fa = sum(response == TRUE)/sum(is_signal == 0), lags = mean(lags))
  calc_creject <- data %>% filter(is_signal == 0) %>% group_by(condition, bin) %>% 
                  dplyr::summarize(creject = sum(response == FALSE)/sum(is_signal == 0), lags = mean(lags))
  calc <- left_join(calc_hit, calc_miss) %>% left_join(calc_fa) %>% left_join(calc_creject)
  calc$hit.correction <- "none"
  calc$fa.correction <- "none"
  calc <- calc %>% select(bin, hit, hit.correction, fa, fa.correction, miss, creject, condition, lags )
  calc
}

test_calc <- satf_calc(data = test_gen_bound)

cmp_dprime <- function(data, hit, fa, miss, creject) {
  n_signal <- hit + miss
  n_noise <- fa + creject
  p_hit <- hit/n_signal
  p_fa <- fa/n_noise
  data$dprimes <- qnorm(p_hit) - qnorm(p_fa)
  data$criterion <- -0.5*(qnorm(p_hit) + qnorm(p_fa))
  drop <- c("miss", "creject", "criterion")
  data <- data[, !(names(data) %in% drop)] %>% 
          select(bin, hit, hit.correction, fa, fa.correction, lags, dprimes, condition) 
  data
}

test_dprimes <- cmp_dprime(data = test_calc, hit = test_calc$hit, miss = test_calc$miss,
                           fa = test_calc$fa, creject = test_calc$creject)
View(test_dprimes)

### their data / start (for comparison)
library(mrsat)

data(Auditory_demo)

my.signal <- list(noint = c(1,3), 
                  obrel = c(5,8), 
                  obrelsub = c(11, 14))
my.noise <- list(noint = c(2,4), 
                 obrel = c(6, 7, 9, 10), 
                 obrelsub = c(12, 13, 15, 16))

d.bins <- get.bins(Auditory_demo, auditory=TRUE)

d.dprime <- get.dprime(d.bins$bins, signal.list=my.signal, noise.list=my.noise, 
                       is.index=TRUE, binmax=14)

View(d.dprime)

### their data / end

pc333 <- list(asym = c(1, 2, 3), rate = c(1, 2, 3), incp = c(1, 2, 3))
fit333 <- fit.SATcurve(test_dprimes, par.cond = pc333)
plot(fit333, main = "333")

pc311 <- list(asym = c(1, 2, 3), rate = c(1, 1, 1), incp = c(1, 1, 1))
fit.311 <- fit.SATcurve(test_dprimes, par.cond = pc311)
plot(fit.311, main="311")

SATsummary.list(list(fit.333, fit.311))