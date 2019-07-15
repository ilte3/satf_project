library(tidyverse)
library(magrittr)
library(ggpubr)

res0_0_0 <- feather::read_feather("../simulations_output/simulation_20-100-0.90-1.50-3.00_0.00-0.00-0.00/all_simulations.feather")
res0_0_0 <- res0_0_0 %>% filter(model == "2-2-2")
View(res0_0_0)

err0_0_0 <- data.frame( err_incp1 = res0_0_0$intercept1 - res0_0_0$true_intercept1,
                        err_incp2 = res0_0_0$intercept2 - res0_0_0$true_intercept2,
                        err_rate1 = res0_0_0$rate1 - res0_0_0$true_rate1,
                        err_rate2 = res0_0_0$rate2 - res0_0_0$true_rate2,
                        err_asymp1 = res0_0_0$asymptote1 - res0_0_0$true_asymptote1,
                        err_asymp2 = res0_0_0$asymptote2 - res0_0_0$true_asymptote2 )

par(mfrow = c(3,2))
hist(err0_0_0$err_incp1, main = "Intercept 1 Error (estimate - true value)")
hist(err0_0_0$err_incp2, main = "Intercept 2 Error (estimate - true value)")
hist(err0_0_0$err_rate1, main = "Rate 1 Error (estimate - true value)")
hist(err0_0_0$err_rate2, main = "Rate 2 Error (estimate - true value)")
hist(err0_0_0$err_asymp1, main = "Asymptote 1 Error (estimate - true value)")
hist(err0_0_0$err_asymp2, main = "Asymptote 2 Error (estimate - true value)")

incp1 <- ggplot(err0_0_0, aes(x = err_incp1)) + geom_histogram()
incp2 <- ggplot(err0_0_0, aes(x = err_incp2)) + geom_histogram()
rate1 <- ggplot(err0_0_0, aes(x = err_rate1)) + geom_histogram()
rate2 <- ggplot(err0_0_0, aes(x = err_rate2)) + geom_histogram()
asymp1 <- ggplot(err0_0_0, aes(x = err_asymp1)) + geom_histogram()
asymp2 <- ggplot(err0_0_0, aes(x = err_asymp2)) + geom_histogram()

fig0_0_0 <- ggarrange(incp1, incp2, rate1, rate2, asymp1, asymp2, ncol = 2, nrow = 3)

fig0_0_0 <- annotate_figure(fig0_0_0, top = text_grob("Error plots - sim0_0_0, model 2-2-2"))
fig0_0_0

res0.05_0_0 <- feather::read_feather("../simulations_output/simulation_20-100-0.90-1.50-3.00_0.05-0.00-0.00/all_simulations.feather")
res0.05_0_0 <- res0.05_0_0 %>% filter(model == "2-2-2")
View(res0.05_0_0)

err0.05_0_0 <- data.frame( err_incp1 = res0.05_0_0$intercept1 - res0.05_0_0$true_intercept1,
                           err_incp2 = res0.05_0_0$intercept2 - res0.05_0_0$true_intercept2,
                           err_rate1 = res0.05_0_0$rate1 - res0.05_0_0$true_rate1,
                           err_rate2 = res0.05_0_0$rate2 - res0.05_0_0$true_rate2,
                           err_asymp1 = res0.05_0_0$asymptote1 - res0.05_0_0$true_asymptote1,
                           err_asymp2 = res0.05_0_0$asymptote2 - res0.05_0_0$true_asymptote2 )

par(mfrow = c(3,2))
hist(err0.05_0_0$err_incp1, main = "Intercept 1 Error (estimate - true value)")
hist(err0.05_0_0$err_incp2, main = "Intercept 2 Error (estimate - true value)")
hist(err0.05_0_0$err_rate1, main = "Rate 1 Error (estimate - true value)")
hist(err0.05_0_0$err_rate2, main = "Rate 2 Error (estimate - true value)")
hist(err0.05_0_0$err_asymp1, main = "Asymptote 1 Error (estimate - true value)")
hist(err0.05_0_0$err_asymp2, main = "Asymptote 2 Error (estimate - true value)")

incp1 <- ggplot(err0.05_0_0, aes(x = err_incp1)) + geom_histogram()
incp2 <- ggplot(err0.05_0_0, aes(x = err_incp2)) + geom_histogram()
rate1 <- ggplot(err0.05_0_0, aes(x = err_rate1)) + geom_histogram()
rate2 <- ggplot(err0.05_0_0, aes(x = err_rate2)) + geom_histogram()
asymp1 <- ggplot(err0.05_0_0, aes(x = err_asymp1)) + geom_histogram()
asymp2 <- ggplot(err0.05_0_0, aes(x = err_asymp2)) + geom_histogram()

fig0.05_0_0 <- ggarrange(incp1, incp2, rate1, rate2, asymp1, asymp2, ncol = 2, nrow = 3)
fig0.05_0_0 <- annotate_figure(fig0.05_0_0, top = text_grob("Error plots - sim0.05_0_0, model 2-2-2"))
fig0.05_0_0

err0.05_0_0 %<>% gather("err_incp1", "err_incp2", value = "err_incp") %>%
                 gather("err_rate1", "err_rate2", value = "err_rate") %>%
                 gather("err_asymp1", "err_asymp2", value = "err_asymp") %<>% 
                 dplyr::select(err_incp, err_rate, err_asymp)

incp_rate_0.05_0_0 <- ggplot(err0.05_0_0, aes(err_incp, err_rate)) + geom_point() + 
  ggtitle("incp_rate") + theme(plot.title = element_text(hjust = 0.5))

incp_asymp_0.05_0_0 <- ggplot(err0.05_0_0, aes(err_incp, err_asymp)) + geom_point() + 
  ggtitle("incp_asymp") + theme(plot.title = element_text(hjust = 0.5))

rate_asymp_0.05_0_0 <- ggplot(err0.05_0_0, aes(err_rate, err_asymp)) + geom_point() + 
  ggtitle("rate_asymp") + theme(plot.title = element_text(hjust = 0.5))

fig0.05_0_0_bivar <- ggarrange(incp_rate_0.05_0_0, incp_asymp_0.05_0_0, rate_asymp_0.05_0_0, ncol = 1, nrow = 3)
fig0.05_0_0_bivar <- annotate_figure(fig0.05_0_0_bivar, 
                                     top = text_grob("Bivariate plots - sim0.05_0_0, model 2-2-2", face = "bold", size = "12"))
fig0.05_0_0_bivar

err0_0_0 %<>% gather("err_incp1", "err_incp2", value = "err_incp") %>%
                 gather("err_rate1", "err_rate2", value = "err_rate") %>%
                 gather("err_asymp1", "err_asymp2", value = "err_asymp") %<>% 
                 dplyr::select(err_incp, err_rate, err_asymp)

incp_rate_0_0_0 <- ggplot(err0_0_0, aes(err_incp, err_rate)) + 
  geom_point(alpha = 0.5) + 
  ggtitle("incp_rate") + 
  theme(plot.title = element_text(hjust = 0.5))

incp_asymp_0_0_0 <- ggplot(err0_0_0, aes(err_incp, err_asymp)) + 
  geom_point(alpha = 0.5) + 
  ggtitle("incp_asymp") + 
  theme(plot.title = element_text(hjust = 0.5))

rate_asymp_0_0_0 <- ggplot(err0_0_0, aes(err_rate, err_asymp)) + 
  geom_point(alpha = 0.5) + 
  ggtitle("rate_asymp") + 
  theme(plot.title = element_text(hjust = 0.5))

fig0_0_0_bivar <- ggarrange(incp_rate_0_0_0, incp_asymp_0_0_0, rate_asymp_0_0_0, ncol = 1, nrow = 3)
fig0_0_0_bivar <- annotate_figure(fig0_0_0_bivar, 
                                     top = text_grob("Bivariate plots - sim0_0_0, model 2-2-2", face = "bold", size = "12"))
fig0_0_0_bivar
