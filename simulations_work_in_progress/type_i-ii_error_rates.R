library(tidyverse)
library(magrittr)

res0_0_0 <- feather::read_feather("../simulations_output/simulation_20-100-0.90-1.50-3.00_0.00-0.00-0.00/all_simulations.feather")
res0_0_0 <- res0_0_0 %>% filter(model == "2-2-2")
View(res0_0_0)

summary(res0_0_0)
