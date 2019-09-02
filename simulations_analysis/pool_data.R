library(tidyverse)
library(magrittr)

process_file <- function(fname)
{
  df <- feather::read_feather(fname)
  
  # reverse model id (the simulation code counterintuitively codes the id as N_asymptotes-N_rates-N_intercepts)
  df$model <- stringi::stri_reverse(df$model)
  
  # extract information 
  sim_id <- fname %>% gsub("/all_simulations.feather$", "", .) %>% gsub("^.*simulation_50-", "", .)
  sim_id_components <- strsplit(sim_id, "[-_]") %>% .[[1]] %>% as.numeric()
  df$sim_n_obs <- sim_id_components[1] %>% as.integer()
  df$sim_avg_intercept <- sim_id_components[2]
  df$sim_avg_rate <- sim_id_components[3]
  df$sim_avg_asymptote <- sim_id_components[4]
  #df$sim_delta_intercept <- sim_id_components[5]
  #df$sim_delta_rate <- sim_id_components[6]
  #df$sim_delta_asymptote <- sim_id_components[7]
  df$delta_intercept <- with(df, true_intercept2 - true_intercept1) %>% round(3)
  df$delta_rate <- with(df, true_rate2 - true_rate1) %>% round(3)
  df$delta_asymptote <- with(df, true_asymptote2 - true_asymptote1) %>% round(3)
  
  df$method <- NULL
  
  df
}

nunique <- function(x) length(unique(x))

select_best_models <- function(df, metric = "adjR2")
{
  df_best <-
  df %>% group_by(#method,
                 simulation_id, participant_id,
                 sim_avg_intercept, sim_avg_rate, sim_avg_asymptote,
                 delta_intercept, delta_rate, delta_asymptote) %>%
    mutate( best_model_adjR2 = (adjR2 == max(adjR2)),
            best_model_AIC = (AIC == min(AIC)) )
  
  if (metric == "adjR2") {
      df_best %<>% filter( best_model_adjR2 )
    
  } else if (metric == "AIC") {
      df_best %<>% filter( best_model_AIC )
    
  } else {
      stop("Unknown metric.")    
  }

  stopifnot( nrow(df_best) == nrow(df %>% dplyr::select(simulation_id, participant_id) %>% unique()) )
  df_best
}


fnames_sr <- dir("../simulations/simulations_output", pattern = ".feather", recursive = T, full.names = T)

df <- plyr::ldply(fnames_sr, function(fname) { process_file(fname) %>% subset(model == "2-2-2") }, .progress = "text")
df$method <- "SR"
feather::write_feather(df, path = "~/sim_sr_222.feather")
gc()

df <- plyr::ldply(fnames_sr, function(fname) { process_file(fname) %>% select_best_models(metric = "adjR2") }, .progress = "text")
df$method <- "SR"
feather::write_feather(df, path = "~/sim_sr_best_adjR2.feather")
gc()

df <- plyr::ldply(fnames_sr, function(fname) { process_file(fname) %>% select_best_models(metric = "AIC") }, .progress = "text")
df$method <- "SR"
feather::write_feather(df, path = "~/sim_sr_best_AIC.feather")
gc()




fnames_mr <- dir("../simulations/simulations_output_mr", pattern = ".feather", recursive = T, full.names = T)

df <- plyr::ldply(fnames_mr, function(fname) { process_file(fname) %>% subset(model == "2-2-2") }, .progress = "text")
df$method <- "MR"
feather::write_feather(df, path = "~/sim_mr_222.feather")
gc()

df <- plyr::ldply(fnames_mr, function(fname) { process_file(fname) %>% select_best_models(metric = "adjR2") }, .progress = "text")
df$method <- "MR"
feather::write_feather(df, path = "~/sim_mr_best_adjR2.feather")
gc()

df <- plyr::ldply(fnames_mr, function(fname) { process_file(fname) %>% select_best_models(metric = "AIC") }, .progress = "text")
df$method <- "MR"
feather::write_feather(df, path = "~/sim_mr_best_AIC.feather")
gc()








df_sr_best <- plyr::ldply(fnames_sr, function(fname) { process_file(fname) %>% subset(model == "2-2-2") }, .progress = "text")
df_sr_222$method <- "SR"
gc()

df_mr <- plyr::ldply(fnames_mr, process_file, .progress = "text")
df_mr$method <- "MR"
gc()

colnames(df_sr)
colnames(df_mr)
df <- rbind(df_sr, df_mr)
feather::write_feather(df, path = "~/simulations_pooled.feather")
gc()

df222 <- subset(df, model == "2-2-2")
feather::write_feather(df222, path = "~/simulations_pooled_222.feather")
gc()


feather::write_feather(df_best, path = "~/simulations_pooled_best.feather")




# head(df)
# tail(df)
# 
# 
# xtabs(~df$sim_avg_intercept)
# xtabs(~df$sim_avg_rate)
# xtabs(~df$sim_avg_asymptote)
# 
# with(df, xtabs(~sim_delta_intercept + delta_intercept))
# with(df, xtabs(~sim_delta_rate + delta_rate))
# with(df, xtabs(~sim_delta_asymptote + delta_asymptote))
