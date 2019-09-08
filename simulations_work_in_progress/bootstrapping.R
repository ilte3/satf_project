library(tidyverse)
library(magrittr)
library(ggpubr)
library(boot)

df <- feather::read_feather(path = "all_simulations.feather")
df222 <- df %>% filter(model == "2-2-2")

fun_t_test <- function(df) {
  df$group_id <- rep(1:(nrow(df)/20), each = 20)
  df %>% group_by(group_id) %>% 
         dplyr::summarize(t_incp = t.test(intercept1, intercept2, paired = T, conf.level = 0.95)$statistic,
                          t_rate = t.test(rate1, rate2, paired = T, conf.level = 0.95)$statistic,
                          t_asymp = t.test(asymptote1, asymptote2, paired = T, conf.level = 0.95)$statistic)
}

t_tests_df <- fun_t_test(df222)

n_obs <- length(t_tests_df)
n_boot <- 1000
t_incp <- t_tests_df$t_incp
t_rate <- t_tests_df$t_rate
t_asymp <- t_tests_df$t_asymp

# bootstrap_incp <- matrix(sample(t_incp, size = n_obs*n_boot, replace = T), nrow = n_obs, ncol = n_boot)
# bootstrap_rate <- matrix(sample(t_rate, size = n_obs*n_boot, replace = T), nrow = n_obs, ncol = n_boot)
# bootstrap_asymp <- matrix(sample(t_asymp, size = n_obs*n_boot, replace = T), nrow = n_obs, ncol = n_boot)

bootstrap_samples <- data.frame(boot_t_incp = sample(t_incp, size = n_obs*n_boot, replace = T),
                                boot_t_rate = sample(t_rate, size = n_obs*n_boot, replace = T),
                                boot_t_asymp = sample(t_asymp, size = n_obs*n_boot, replace = T))

###########################################################################################################