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
  # stopifnot(nrow(df)/20 == round(nrow(df)/20)) 
  df <- df %T>% {.$n_for_t_tests = rep(c(1:(nrow(dfff)/20)), each = 20)}
}
dfff <- conv_n_participants(dfff)

#The first 20 would have an output like this 1-row data frame below (t_tests_df).

t_test_incp <- t.test(dfff$intercept1[0:20], dfff$intercept2[0:20], paired = T, conf.level = 0.95)
t_test_rate <- t.test(dfff$rate1[0:20], dfff$rate2[0:20], paired = T, conf.level = 0.95)
t_test_asymp <- t.test(dfff$asymptote1[0:20], dfff$asymptote2[0:20], paired = T, conf.level = 0.95)


t_tests_df <- data.frame(t_incp = t_test_incp$statistic[["t"]], t_rate = t_test_rate$statistic[["t"]],
                         t_asymp = t_test_asymp$statistic[["t"]])

t_tests_df <- t_tests_df %T>% {.$dif_significance_incp = ifelse(.$t_incp > 2, "YES", "NO")} %T>%
                              {.$dif_significance_rate = ifelse(.$t_rate > 2, "YES", "NO")} %T>%
                              {.$dif_significance_asymp = ifelse(.$t_asymp > 2, "YES", "NO")} %>%
                              dplyr::select(t_incp, dif_significance_incp, t_rate, dif_significance_rate,
                                           t_asymp, dif_significance_asymp)
t_tests_df

xxx <- plyr::ldply(1:dfff$n_for_t_tests, function(df) {
    t_test_incp <- t.test(dfff$intercept1, dfff$intercept2, paired = T, conf.level = 0.95)
    t_test_rate <- t.test(dfff$rate1, dfff$rate2, paired = T, conf.level = 0.95)
    t_test_asymp <- t.test(dfff$asymptote1, dfff$asymptote2, paired = T, conf.level = 0.95)
    
    t_tests_df <- data.frame(t_incp = t_test_incp$statistic[["t"]], t_rate = t_test_rate$statistic[["t"]],
                             t_asymp = t_test_asymp$statistic[["t"]])
    
    t_tests_df <- t_tests_df %T>% {.$dif_significance_incp = ifelse(.$t_incp > 2, "YES", "NO")} %T>%
      {.$dif_significance_rate = ifelse(.$t_rate > 2, "YES", "NO")} %T>%
      {.$dif_significance_asymp = ifelse(.$t_asymp > 2, "YES", "NO")} %>%
      dplyr::select(t_incp, dif_significance_incp, t_rate, dif_significance_rate,
                    t_asymp, dif_significance_asymp)
    t_tests_df
})

for (i in dfff$n_for_t_tests) {
  t_test_incp <- t.test(dfff$intercept1[i], dfff$intercept2[i], paired = T, conf.level = 0.95)
  t_test_rate <- t.test(dfff$rate1[i], dfff$rate2[i], paired = T, conf.level = 0.95)
  t_test_asymp <- t.test(dfff$asymptote1[i], dfff$asymptote2[i], paired = T, conf.level = 0.95)
}


