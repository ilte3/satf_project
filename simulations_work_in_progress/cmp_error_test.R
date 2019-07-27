library(tidyverse)
library(magrittr)
library(ggpubr)

res0_0_0 <- feather::read_feather("../simulations/simulations_output/simulation_20-100-0.90-1.50-3.00_0.00-0.00-0.00/all_simulations.feather")
res0_0_0 <- res0_0_0 %>% filter(model == "2-2-2")
head(res0_0_0)

est <- res0_0_0 #%>% dplyr::select(-R2, -adjR2, -logLik, -AIC)
est1 <- est %>% dplyr::rename(#simulation_id, model, participant_id,
                              true_intercept=true_intercept1, true_rate=true_rate1, true_asymptote=true_asymptote1, 
                              est_intercept=intercept1, est_rate=rate1, est_asymptote=asymptote1) %>% 
                dplyr::select(-true_intercept2, -true_rate2, -true_asymptote2, 
                              -intercept2, -rate2, -asymptote2)
est2 <- est %>% dplyr::rename(#simulation_id, model, participant_id,
                              true_intercept=true_intercept2, true_rate=true_rate2, true_asymptote=true_asymptote2, 
                              est_intercept=intercept2, est_rate=rate2, est_asymptote=asymptote2) %>% 
                dplyr::select(-true_intercept1, -true_rate1, -true_asymptote1, 
                              -intercept1, -rate1, -asymptote1)
est1$condition <- 1
est2$condition <- 2
est <- dplyr::bind_rows(est1, est2) %>% dplyr::arrange(simulation_id, model, participant_id, condition)

head(est) %>% as.data.frame()

  
est %<>% dplyr::mutate(err_intercept = est_intercept - true_intercept,
                       err_rate = est_rate - true_rate,
                       err_asymp = est_asymptote - true_asymptote)

head(est) %>% as.data.frame()

ggplot(est, aes(err_intercept)) + geom_histogram()
ggplot(est, aes(err_rate)) + geom_histogram()
ggplot(est, aes(err_asymp)) + geom_histogram()

# check out the bump in the histogram below -0.2
est %>% filter(round(err_intercept,2) == -0.2) %>% summary()

est %>% subset(err_intercept < -0.1 & err_intercept > -0.4) %>% 
  ggplot(aes(err_intercept, err_rate)) + 
  geom_point(alpha = 0.05)

windows()

ggplot(est, aes(err_intercept, err_rate)) + 
      geom_point(alpha = 0.05)

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



###
### check which model is selected most often
###

res <- feather::read_feather("../simulations/simulations_output/simulation_20-100-0.90-1.50-3.00_0.10-0.00-0.00/all_simulations.feather")

head(res)


res$diff_intercept <- with(res, grepl("^2", model) )
res$diff_rate <- with(res, grepl("^.-2", model) )
res$diff_asymptote <- with(res, grepl("^.-.-2", model) )

x <- 
res %>% group_by(simulation_id, participant_id) %>% 
        dplyr::mutate(best = adjR2 == max(adjR2)) %>%
        filter(best) %>% group_by(model, diff_intercept, diff_rate, diff_asymptote) %>%
        dplyr::summarise(N = length(model))
x$P <- 100*x$N/sum(x$N)
x

x %>% subset(diff_intercept) %>% .$P %>% sum()
x %>% subset(diff_rate) %>% .$P %>% sum()
x %>% subset(diff_intercept | diff_rate) %>% .$P %>% sum()
x %>% subset(diff_asymptote) %>% .$P %>% sum()

x %>% group_by(diff_dynamics = diff_intercept | diff_rate, diff_asymptote) %>% dplyr::summarise(P = sum(P))

res0_0_0 <- res0_0_0 %>% filter(model == "2-2-2")
head(res0_0_0)

