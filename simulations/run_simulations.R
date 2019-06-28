#source("./simulations_functions.R")
library(mrsat)
#library(doMC)
#doMC::registerDoMC(4)

run_simulation <- function(n_simulations, n_participants, n_trials_per_interval, avg_incp, avg_rate, avg_asymp, delta_intercept, delta_rate, delta_asymptote) 
{
  dir <- sprintf("./simulations_output/simulation_%d-%d-%0.2f-%0.2f-%0.2f_%0.2f-%0.2f-%0.2f", n_participants, n_trials_per_interval, avg_incp, avg_rate, avg_asymp,  delta_intercept, delta_rate, delta_asymptote)
  if (!file.exists(dir))
    dir.create(dir, recursive = T)
  
  intercept <- avg_incp + c(-.5, .5)*delta_intercept
  rate <- avg_rate + c(-.5, .5)*delta_rate
  asymptote <- avg_asymp + c(-.5, .5)*delta_asymptote
  
  intercepts <- rep(list(intercept), n_participants)
  rates <- rep(list(rate), n_participants)
  asymptotes <- rep(list(asymptote), n_participants)
  
  sim_participants <- function(n_trials_per_interval, time, intercepts, rates, asymptotes, simulation_id)
  {
      estimates <- 
      plyr::ldply(1:n_participants, function(i) {
          estimates_i <- sim_participant(n_trials_per_interval, time, intercepts[[i]], rates[[i]], asymptotes[[i]], 
                                         debug_fname = paste0(dir, "/", simulation_id, "_", i, ".rda") )
          estimates_i %T>% {.$participant_id = i}
      }, .progress = "text")
      
      estimates <- estimates %T>% {.$simulation_id <- simulation_id} %>% 
        dplyr::select(simulation_id, model, participant_id, 
                      true_intercept1, true_intercept2, 
                      true_rate1, true_rate2, 
                      true_asymptote1, true_asymptote2, 
                      R2, adjR2, logLik, AIC, 
                      intercept1, intercept2, 
                      rate1, rate2, 
                      asymptote1, asymptote2)
      
      feather::write_feather( estimates, path = paste0(dir, "/", simulation_id, ".feather") )

      estimates
  }
  
  sim_res <- plyr::ldply(1:n_simulations, function(j) {
      estimates <- sim_participants(n_trials_per_interval, time, intercepts, rates, asymptotes, j)
  }, .parallel = F)# , .progress = "text")
  
  
  feather::write_feather(sim_res, path = paste0(dir, "/", "all_simulations", ".feather"))
}


###################################################################################################################################

#options(warn = 1)

n_simulations <- 1000
n_participants <- 20
n_trials_per_interval <- 100

time <- seq(0, 5.6, 0.35)

avg_incp <- .9
avg_rate <- 1.5
avg_asymp <- 3.0

run_simulation(n_simulations, n_participants, n_trials_per_interval, avg_incp, avg_rate, avg_asymp, 
               delta_intercept = 0.0, delta_rate = 0.0, delta_asymptote = 0.0) 
run_simulation(n_simulations, n_participants, n_trials_per_interval, avg_incp, avg_rate, avg_asymp, 
               delta_intercept = 0.05, delta_rate = 0.0, delta_asymptote = 0.0) 
run_simulation(n_simulations, n_participants, n_trials_per_interval, avg_incp, avg_rate, avg_asymp, 
               delta_intercept = 0.1, delta_rate = 0.0, delta_asymptote = 0.0) 

res0_0_0 <- feather::read_feather("./simulations_output/simulation_20-100-0.90-1.50-3.00_0.00-0.00-0.00/all_simulations.feather") 
res0.05_0_0 <- feather::read_feather("./simulations_output/simulation_20-100-0.90-1.50-3.00_0.05-0.00-0.00/all_simulations.feather") 
res0.10_00 <- feather::read_feather("./simulations_output/simulation_20-100-0.90-1.50-3.00_0.10-0.00-0.00/all_simulations.feather") 


summary(res0_0_0)
head(res0_0_0)


subset(res0_0_0, is.na(R2)) %>% View()
