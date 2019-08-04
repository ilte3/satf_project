library(tidyverse)
library(magrittr)
library(ggpubr)

# Using any all_simulations.feather from the output folder.
# I used a pooled version for other analyses but
# since this is trial and error, a smaller df would suffice.

df <- feather::read_feather(path = "C:/Users/Ýlte/Desktop/all_simulations.feather")
df222 <- df %>% filter(model == "2-2-2")

##### err

# err1 <- with(df,data.frame(err_incp = intercept1 - true_intercept1,
#                            err_rate = rate1 - true_rate1,
#                            err_asymp = asymptote1 - true_asymptote1))
# 
# err2 <- with(df,data.frame(err_incp = intercept2 - true_intercept2,
#                            err_rate = rate2 - true_rate2,
#                            err_asymp = asymptote2 - true_asymptote2))
# 
# err <- rbind(err1, err2)
#
# incp_rate <- ggplot(err, aes(err_incp, err_rate)) +
#   geom_point(alpha = 0.5) +
#   ggtitle("incp_rate") +
#   theme(plot.title = element_text(hjust = 0.5))
#
# incp_asymp <- ggplot(err, aes(err_incp, err_asymp)) +
#   geom_point(alpha = 0.5) +
#   ggtitle("incp_asymp") +
#   theme(plot.title = element_text(hjust = 0.5))
#
# rate_asymp <- ggplot(err, aes(err_rate, err_asymp)) +
#   geom_point(alpha = 0.5) +
#   ggtitle("rate_asymp") +
#   theme(plot.title = element_text(hjust = 0.5))
#
# fig_bivar <- ggarrange(incp_rate, incp_asymp, rate_asymp, ncol = 1, nrow = 3)
# fig_bivar <- annotate_figure(fig_bivar,
#                              top = text_grob("Bivariate plots - mode: SR - model 2-2-2", face = "bold", size = "12"))
# fig_bivar



##### t-test

#I would like to create a loop that repeats the t-test function every 20 row.
dfff <- df222
conv_n_participants <- function(df) {
  df <- df %T>% {.$n_participants_experimental = rep(c(1:20), nrow(df)/20)}
}
dfff <- conv_n_participants(dfff)

#The first 20 would have an output like this 1-row data frame below.
t_tests <- t.test(df222$intercept1[0:20], df222$intercept2[0:20], paired = T, conf.level = 0.95)
library(tidystats)
t_tests_df <- tidy_stats(t_tests)
t_tests_df <- with(t_tests_df, data.frame(t_value = value[2], p_value = value[4], ci_lower = value[5], ci_upper = value[6]))
t_tests_df$dif_significance <- with(t_tests_df, ifelse(t_value > 2, "YES", "NO"))
