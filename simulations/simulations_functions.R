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


satf_gen_cond_mr <- function(mu, criterion, time, n, label) {
  stopifnot(length(mu) == length(criterion) && length(criterion) == length(time))
  trial = rep(1:n, each = length(time))
  interval <- 1:length(time)
  
  df <- data.frame(condition = label, interval = interval, time = time, trial = trial, mu = mu, 
                   criterion = criterion, delta_evidence = c(0, diff(mu)) )
  
  # generate start position (for index=1), and delta due to noise (for index > 1) 
  df$evidence_rel_position <- rnorm(nrow(df), mean = 0, sd = 1)
  # determine the relative evidence positions (modulo the provided mus)
  df %<>% group_by(trial) %>% dplyr::mutate(evidence_rel_position = cumsum(evidence_rel_position)/sqrt(interval),
                                            evidence_position = evidence_rel_position + mu#,
                                            #cor = cor(evidence_rel_position, lag(evidence_rel_position), use = "complete.obs")
  )
  # NOTE: The above assumptions about the generative process for MR SAT data results in the serial correlation between
  #       evidence_rel_position values being abound ~0.6. Other assumptions (lower or higher SD) would result in a different
  #       correlation coefficient.
  
  df$response <- df$evidence_position > criterion
  df %>% ungroup() %>% dplyr::select(condition, interval, time, trial, response)
}

satf_gen_cond <- function(mu, criterion, time, n, label) {
  stopifnot(length(mu) == length(criterion) && length(criterion) == length(time))
  trial = rep(1:n, each = length(time))
  interval <- 1:length(time)
  df <- data.frame(condition = label, interval = interval, time = time, trial = trial, mu = mu, criterion = criterion)
  df$evidence_position <- rnorm(nrow(df), mean = mu, sd = 1)
  df$response <- df$evidence_position > criterion
  df %>% dplyr::select(condition, interval, time, trial, response)
}

satf_gen <- function(time, n, intercept, rate, asymptote, label = "condition1", fn_satf_gen_cond = satf_gen_cond) {
  dprime <- satf(time, intercept, rate, asymptote)
  criterion <- 0.5*dprime
  df0 <- fn_satf_gen_cond(mu = dprime*0, criterion = criterion, time = time, n = n, label = label)
  df1 <- fn_satf_gen_cond(mu = dprime, criterion = criterion, time = time, n = n, label = label)
  df0$is_signal <- 0
  df1$is_signal <- 1
  df1$trial = df1$trial + max(df0$trial)
  rbind(df0, df1)
}

satf_gen_mr <- function(time, n, intercept, rate, asymptote, label) {
  satf_gen(time = time, n = n, intercept = intercept, rate = rate, asymptote = asymptote, label = label, fn_satf_gen_cond = satf_gen_cond_mr)
}

cmp_acc_stat <- function(data)
{
  calc <- data %>% group_by(condition, interval, is_signal) %>% 
                   dplyr::summarize(p_yes = mean(response), time = mean(time), n = length(response))
  calc$correction <- 'none'

  # 'correct' extreme values
  is_zero <- calc$p_yes == 0
  is_one <- calc$p_yes == 1
  
  if (any(is_zero) || any(is_one)) {
    calc$p_yes[is_zero] <- with(calc[is_zero,], 0.5/n)
    calc$p_yes[is_one] <- with(calc[is_one,], (n-0.5)/n)
    calc$correction[is_zero | is_one] <- 'extreme'
  }
  
  calc_signal <- calc %>% filter(is_signal == 1) %>% 
                    dplyr::select(condition, interval, time, hit = p_yes, 
                                  n_signal = n, hit.correction = correction) %>% 
                    dplyr::mutate(miss = 1-hit)
  calc_noise <- calc %>% filter(is_signal == 0) %>% 
                    dplyr::select(condition, interval, time, fa = p_yes, 
                                  n_noise = n, fa.correction = correction) %>% 
                    dplyr::mutate(creject = 1-fa)

  left_join(calc_signal, calc_noise, by = c("condition", "interval", "time")) %>% 
            dplyr::select(interval, time, hit, miss, hit.correction, n_signal, 
                                          fa, creject, fa.correction, n_noise, 
                          condition )
}

cmp_dprime <- function(data) {
  n_signal <- data$hit + data$miss
  n_noise <- data$fa + data$creject
  p_hit <- data$hit/n_signal
  p_fa <- data$fa/n_noise
  
  data$dprime <- qnorm(p_hit) - qnorm(p_fa)
  data$criterion <- -0.5*(qnorm(p_hit) + qnorm(p_fa))
  data %<>% dplyr::select(interval, time, hit, miss, hit.correction, n_signal, 
                                          fa, creject, fa.correction, n_noise, 
                          dprime, criterion, condition) 
  data
}

mrsat_fitcurve <- function(data, pc = list(asym = c(1, 2), rate = c(1, 2), incp = c(1, 2)), 
                                 startvals, show_plot = FALSE)
{
  data <- data %>%  dplyr::select(bin = interval, hit, hit.correction,
                                  hit.denom = n_signal, fa, fa.correction, fa.denom = n_noise, 
                                  lags = time, dprimes = dprime, condition)
  
  par_constraints = get.param(pc, auto.asym = FALSE, 
                              asym = c(startvals$asymptote, 0.1, 5), 
                              rate = c(startvals$rate, 0.3, 5), 
                              incp = c(startvals$intercept, 0.1, 2), data = data)
  fit <- fit.SATcurve(data, par.cond = pc, params = par_constraints, maxit = 10^6)
  sum_curve <- summary.SATcurve(fit)
  convergence <- fit$fit[c('convergence', 'iterations')]
  sum_curve %<>% cbind( convergence ) 
  
  n_unique <- function(x) length(unique(x))
  
  model_id <- sapply(pc, n_unique) %>% paste(collapse = "-")

  model_fit <- sum_curve#[c("incp1", "incp2", "rate1", "rate2", "asym1", "asym2", "R2", "adjR2", "logLik", "AIC")]
  model_fit$model <- model_id
  model_fit %<>% dplyr::rename(intercept1 = incp1, intercept2 = incp2, asymptote1 = asym1, asymptote2 = asym2)
  
  model_fit
}

sim_participant_mr <- function(n, time, intercept, rate, asymptote, fit_start, debug_fname = NULL) {
  sim_participant(n, time, intercept, rate, asymptote, fit_start, debug_fname = NULL, fn_satf_gen = satf_gen_mr) 
    
}

sim_participant <- function(n, time, intercept, rate, asymptote, fit_start,
                            debug_fname = NULL, fn_satf_gen = satf_gen) 
{
  responses1 <- fn_satf_gen(time = time, n = n, intercept = intercept[1], rate = rate[1], asymptote = asymptote[1], label = "condition1")
  responses2 <- fn_satf_gen(time = time, n = n, intercept = intercept[2], rate = rate[2], asymptote = asymptote[2], label = "condition2")
  responses <- rbind(responses1, responses2)
  
  acc_stat <- cmp_acc_stat(data = responses)
  acc_stat <- cmp_dprime(data = acc_stat)
  
  res111 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 1), rate = c(1, 1), incp = c(1, 1)), fit_start)
  res112 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 1), rate = c(1, 1), incp = c(1, 2)), fit_start)
  res121 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 1), rate = c(1, 2), incp = c(1, 1)), fit_start)
  res122 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 1), rate = c(1, 2), incp = c(1, 2)), fit_start)
  res211 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 2), rate = c(1, 1), incp = c(1, 1)), fit_start)
  res221 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 2), rate = c(1, 2), incp = c(1, 1)), fit_start)
  res212 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 2), rate = c(1, 1), incp = c(1, 2)), fit_start)
  res222 <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 2), rate = c(1, 2), incp = c(1, 2)), fit_start)
  res <- dplyr::bind_rows(res111, res112, res121, res122, res211, res221, res212, res222)
  
  true_params <- data.frame(true_intercept1 = intercept[1], true_intercept2 = intercept[2],
                            true_rate1 = rate[1], true_rate2 = rate[2],
                            true_asymptote1 = asymptote[1], true_asymptote2 = asymptote[2])
  res %<>% cbind(true_params)
  res %<>% dplyr::select(model, true_intercept1:true_asymptote2, asymptote1:intercept2, R2:iterations)
  
  attr(res, "responses") <- responses  
  attr(res, "acc_stat") <- acc_stat  
  
  if ( !is.null(debug_fname) && any(is.na(res$R2)) ) {
      save(responses, res, file = debug_fname)
  }
   
  res
}
