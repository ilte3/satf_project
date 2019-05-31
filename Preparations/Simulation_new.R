source("../default.R")
library(mrsat)

###################################################################################################################################

options(warn = 1)

n_simulations <- 1
n_participants <- 1
n_trials_per_interval <- 50

time <- seq(0, 5.6, 0.35)

avg_incp <- rnorm(1, 0.4, 0.1)^2
avg_rate <- rnorm(1, 0.4, 0.1)^2
avg_asymp <- rnorm(1, 3, 1)

delta_intercept = 0.05
delta_rate = 0.05
delta_asymptote = 0.05

intercept <- avg_incp + c(-.5, .5)*delta_intercept
rate <- avg_rate + c(-.5, .5)*delta_rate
asymptote <- avg_asymp + c(-.5, .5)*delta_asymptote

intercepts <- rep(list(intercept), n_participants)
rates <- rep(list(rate), n_participants)
asymptotes <- rep(list(asymptote), n_participants)

sim_tmp <- plyr::ldply(1:n_simulations, function(j) {
  estimates_j <- plyr::ldply(1:n_participants, function(i) {
    estimates_i <- sim_participant(n_trials_per_interval, time,
                                 intercepts[[i]], rates[[i]], asymptotes[[i]], 
                                 show_plot = FALSE)
    estimates_i %T>% {.$participant_id = i}
  })
  estimates_j %T>% {.$simulation_id <- j} %>% 
    dplyr::select(simulation_id, model, participant_id, true_intercept1, true_intercept2, true_rate1, true_rate2, 
                  true_asymptote1, true_asymptote2, R2, adjR2, AIC, intercept1, intercept2, rate1, rate2, asymptote1, asymptote2)
})

