
library(tidyverse)
library(magrittr)

#df <- feather::read_feather(path = "~/simulations_pooled.feather")

df_best <- feather::read_feather(path = "~/simulations_pooled_best.feather")


df_summary <-
  df_best %>% group_by(mode, 
                     sim_avg_intercept, sim_avg_rate, sim_avg_asymptote,
                     delta_intercept, delta_rate, delta_asymptote,
                     model) %>%
            summarize(N = n(), 
                      sign_delta_intercept = mean((intercept2-intercept1) > 0),
                      sign_delta_rate = mean((rate2-rate1) > 0),
                      sign_delta_asymptote = mean((asymptote2-asymptote1) > 0) ) %>% 
            mutate(P = N / sum(N))


p_2params <-
df_summary %>% group_by(mode,
                        sim_avg_intercept, sim_avg_rate, sim_avg_asymptote,
                        delta_intercept, delta_rate, delta_asymptote) %>%
                mutate(model_2icpt = grepl("^2", model), 
                       model_2rate = grepl("^[12]-2", model), 
                       model_2asym = grepl("2$", model)
                      ) %>%
                summarize(p_2icpt = sum( P[model_2icpt] ),
                          sgnpos_2icpt = sum( P[model_2icpt]/sum(P[model_2icpt]) * sign_delta_intercept[model_2icpt] ),

                          p_2rate = sum( P[model_2rate]),
                          sgnpos_2rate = sum( P[model_2rate]/sum(P[model_2rate]) * sign_delta_rate[model_2rate] ),
                          
                          p_2asym = sum( P[model_2asym]),
                          sgnpos_2asym = sum( P[model_2asym]/sum(P[model_2asym]) * sign_delta_asymptote[model_2asym] )
                         )

head(p_2params, 20) %>% as.data.frame() 





library(ggplot2)

p_2params %>% ggplot(aes(delta_intercept, `P[2 icpt]`, color = as.factor(delta_rate))) + 
              geom_point(aes(shape = mode)) + geom_line(aes(group = paste(delta_rate,mode))) +
              facet_grid(sim_avg_intercept~delta_asymptote)


View(df_mr)

df_best_mr <- subset(df_best, mode == "MR")  
df_best_mr

df_mr <- subset(df, mode == "MR")
df_mr222 <- subset(df_mr, model == "2-2-2")

df_mr222 %<>% mutate(error_asymptote1 = asymptote1-true_asymptote1, 
                     error_asymptote2 = asymptote2-true_asymptote2, 
                     error_rate1 = rate1-true_rate1, 
                     error_rate2 = rate2-true_rate2, 
                     error_intercept1 = intercept1-true_intercept1, 
                     error_intercept2 = intercept2-true_intercept2,
                     
                     error_delta_asymptote = (true_asymptote2 - true_asymptote1) - (asymptote2 - asymptote1),
                     error_delta_rate = (true_rate2 - true_rate1) - (rate2 - rate1),
                     error_delta_intercept = (true_intercept2 - true_intercept1) - (intercept2 - intercept1)
                    )

df_mr222 %>% ggplot(aes(error_delta_intercept)) + geom_histogram()


time = seq(0,5, .5)
xsr <- satf_gen(time, n=60, intercept=.8, rate=1, asymptote=2, label = "xxx")
xmr <- satf_gen_mr(time, n=60, intercept=.8, rate=1, asymptote=2, label = "xxx")

xsr %>% cmp_acc_stat() %>% cmp_dprime()
xmr %>% cmp_acc_stat() %>% cmp_dprime()

head(xsr)
head(xmr)
