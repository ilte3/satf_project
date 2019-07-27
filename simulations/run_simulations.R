source("./simulations_functions.R")
library(mrsat)
library(doMC)
doMC::registerDoMC(parallel::detectCores(logical=FALSE))

simulate_mrsat = TRUE

run_simulation <- function(n_simulations, n_participants, n_trials_per_interval, avg_incp, avg_rate, avg_asymp, delta_intercept, delta_rate, delta_asymptote) 
{
  if (simulate_mrsat)
    fname_template <- "./simulations_output_mr/simulation_%d-%d-%0.3f-%0.3f-%0.3f_%0.3f-%0.3f-%0.3f"
  else
    fname_template <- "./simulations_output/simulation_%d-%d-%0.3f-%0.3f-%0.3f_%0.3f-%0.3f-%0.3f"
  
  dir <- sprintf(fname_template, 
                 n_participants, n_trials_per_interval, 
                 avg_incp, avg_rate, avg_asymp, 
                 delta_intercept, delta_rate, delta_asymptote)
  fname_all <- paste0(dir, "/", "all_simulations", ".feather")
  
  if (file.exists(dir)) {
      return (NULL)
  } else {
      dir.create(dir, recursive = T)
  }
 
  intercept <- avg_incp + c(-.5, .5)*delta_intercept
  rate <- avg_rate + c(-.5, .5)*delta_rate
  asymptote <- avg_asymp + c(-.5, .5)*delta_asymptote
  
  intercepts <- rep(list(intercept), n_participants)
  rates <- rep(list(rate), n_participants)
  asymptotes <- rep(list(asymptote), n_participants)
  
  if (simulate_mrsat)
    fn_sim_participant <- sim_participant_mr
  else
    fn_sim_participant <- sim_participant

  sim_participants <- function(n_trials_per_interval, time, intercepts, rates, asymptotes, simulation_id)
  {
      estimates <- 
      plyr::ldply(1:n_participants, function(i) {
          estimates_i <- fn_sim_participant(n_trials_per_interval, time, intercepts[[i]], rates[[i]], asymptotes[[i]], 
                                            fit_start = list(intercept = avg_incp, rate = avg_rate, asymptote = avg_asymp),
                                            debug_fname = paste0(dir, "/", simulation_id, "_", i, ".rda") )
          estimates_i %T>% {.$participant_id = i}
      }) #, .progress = "text")

      estimates <- estimates %T>% {.$simulation_id <- simulation_id} %>% 
                    dplyr::select(simulation_id, model, participant_id,
                                  true_intercept1:iterations)

      save( estimates, file = paste0(dir, "/", simulation_id, ".rda") )

      estimates
  }
  
  sim_res <- plyr::ldply(1:n_simulations, function(j) {
      estimates <- sim_participants(n_trials_per_interval, time, intercepts, rates, asymptotes, j)
  }, .parallel = T) #, .progress = "text")
  
  feather::write_feather(sim_res, path = fname_all)
}


###################################################################################################################################

#options(warn = 1)

n_simulations <- 100
n_participants <- 50
time <- seq(0, 5.6, 0.35)

run_cur_sim <- function(n_obs_per_interval_per_cond, avg_incp, avg_rate, avg_asymp, delta_intercept, delta_rate, delta_asymptote) {
    run_simulation(n_simulations, n_participants, n_obs_per_interval_per_cond, avg_incp, avg_rate, avg_asymp, 
                   delta_intercept = delta_intercept, delta_rate = delta_rate, delta_asymptote = delta_asymptote) 
}

for (n_obs_per_interval_per_cond in c(80, 60, 40, 20)) {
  for (avg_rate in c(1, 1.2, 0.8) ) {
    for (avg_asymp in c(2, 1, 3)  ) {
      for (avg_incp in (.7 + .350*c(3/4, 2/4, 1/4, 0)) ) {
        for (delta_asymptote in c(0, .125, .25, .5, 1)) {
          for (delta_rate in c(0, .025, .05, .1, .2)) {
            for (delta_intercept in c(0, .025, .05, .1, .15)) {
              run_cur_sim(n_obs_per_interval_per_cond, avg_incp, avg_rate, avg_asymp, delta_intercept, delta_rate, delta_asymptote)
            }
          }
        }
      }
    }
  }
}
