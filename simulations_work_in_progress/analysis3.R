library(tidyverse)
library(magrittr)
library(ggpubr)

df <- feather::read_feather(path = "all_simulations.feather")
#df <- feather::read_feather(path = "../simulations/simulations_output_mr/simulation_50-20-0.700-0.800-3.000_0.000-0.000-0.000/all_simulations.feather")
df222 <- df %>% filter(model == "2-2-2")

n_participant_exp <- function(df) {
  group <- gl((nrow(df)/20), 20)
  spl_df <- split(df, group)
}

fun_t_test <- function(df) {
  df <- with(df, data.frame(t_incp = t.test(intercept1, intercept2, paired = T, conf.level = 0.95)$statistic,
                            t_rate = t.test(rate1, rate2, paired = T, conf.level = 0.95)$statistic,
                            t_asymp = t.test(asymptote1, asymptote2, paired = T, conf.level = 0.95)$statistic))
  
  df <- df %T>% {.$dif_significance_incp = ifelse(.$t_incp > 2, "YES", "NO")} %T>%
                {.$dif_significance_rate = ifelse(.$t_rate > 2, "YES", "NO")} %T>%
                {.$dif_significance_asymp = ifelse(.$t_asymp > 2, "YES", "NO")} %>%
                dplyr::select(t_incp, dif_significance_incp, t_rate, dif_significance_rate,
                              t_asymp, dif_significance_asymp)
}

################################################################################################################

df222 <- n_participant_exp(df222)
t_tests_df <- as.data.frame(t(sapply(df222, FUN = fun_t_test)))


################################################################################################################
########### And this is what it should have been 
################################################################################################################

df222$group_id <- rep(1:(nrow(df222)/20), each = 20)

t_tests_df <- 
df222 %>% group_by(group_id) %>% 
          dplyr::summarize( t_incp = t.test(intercept1, intercept2, paired = T, conf.level = 0.95)$statistic,
                            t_rate = t.test(rate1, rate2, paired = T, conf.level = 0.95)$statistic,
                            t_asymp = t.test(asymptote1, asymptote2, paired = T, conf.level = 0.95)$statistic )
