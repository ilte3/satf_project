res0_0_0 <- feather::read_feather("../simulations_output/simulation_20-100-0.90-1.50-3.00_0.00-0.00-0.00/all_simulations.feather")
res0_0_0 <- res0_0_0 %>% filter(model == "2-2-2")
View(res0_0_0)

res0_0.05_0 <- feather::read_feather("../simulations_output/simulation_20-100-0.90-1.50-3.00_0.00-0.05-0.00/all_simulations.feather")
res0_0.05_0 <- res0_0.05_0 %>% filter(model == "2-2-2")
View(res0_0.05_0)

res0.05_0_0 <- feather::read_feather("../simulations_output/simulation_20-100-0.90-1.50-3.00_0.05-0.00-0.00/all_simulations.feather")
res0.05_0_0 <- res0.05_0_0 %>% filter(model == "2-2-2")
View(res0.05_0_0)

res0.05_0.05_0 <- feather::read_feather("../simulations_output/simulation_20-100-0.90-1.50-3.00_0.05-0.05-0.00/all_simulations.feather")
res0.05_0.05_0 <- res0.05_0.05_0 %>% filter(model == "2-2-2")
View(res0.05_0.05_0)

res0.10_0_0 <- feather::read_feather("../simulations_output/simulation_20-100-0.90-1.50-3.00_0.10-0.00-0.00/all_simulations.feather")
res0.10_0_0 <- res0.10_0_0 %>% filter(model == "2-2-2")
View(res0.10_0_0)

#res0.10_0.05_0 <- feather::read_feather("../simulations_output/simulation_20-100-0.90-1.50-3.00_0.10-0.05-0.00/all_simulations.feather")

res0.15_0_0 <- feather::read_feather("../simulations_output/simulation_20-100-0.90-1.50-3.00_0.15-0.00-0.00/all_simulations.feather")
res0.15_0_0 <- res0.15_0_0 %>% filter(model == "2-2-2")
View(res0.15_0_0)

res0.20_0_0 <- feather::read_feather("../simulations_output/simulations_output/simulation_20-100-0.90-1.50-3.00_0.20-0.00-0.00/all_simulations.feather")
res0.20_0_0 <- res0.20_0_0 %>% filter(model == "2-2-2")
View(res0.20_0_0)

res0.30_0_0 <- feather::read_feather("../simulations_output/simulations_output/simulation_20-100-0.90-1.50-3.00_0.30-0.00-0.00/all_simulations.feather")
res0.30_0_0 <- res0.30_0_0 %>% filter(model == "2-2-2")
View(res0.30_0_0)

summary(res0_0_0)
summary(res0_0.05_0)
summary(res0.05_0_0)
summary(res0.05_0.05_0)
summary(res0.10_0_0)
summary(res0.15_0_0)
summary(res0.20_0_0)
summary(res0.30_0_0)