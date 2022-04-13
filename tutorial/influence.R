# Install needed packages
install.packages("tidyverse")
install.packages("lme4")
install.packages("lmerTest")

# R Code: Load data
library(tidyverse)
edge = read_csv("./edge_influence.csv")
node = read_csv("./node_influence.csv")


# R Code: Join node data to edge by source id
full_edge_data = left_join(edge,node, by = c("source_ID" = "node_ID"))


# R Code: Calculate Dyadic exposure
full_edge_data = full_edge_data %>% mutate(exposure = communicaiton_count * negative_affect_0)

# Calculate mean exposure 
mean_exposure_tbl = full_edge_data %>% group_by(target_ID) %>% summarize(exposure_mean = mean(exposure,na.rm = TRUE))

# Include mean exposure in the node's data and handle cases where there was no exposure
node = left_join(node,mean_exposure_tbl, by = c("node_ID" = "target_ID"))
node = node %>% mutate(exposure_mean = replace_na(exposure_mean, 0))


# R Code: Run OLS analysis
model_OLS = lm(negative_affect_1 ~ negative_affect_0 + exposure_mean, data = node)

# R Code: Run Multi-Level Analysis accounting for random intercepts and random slopes
library(lme4)
library(lmerTest)
model_int = lmer(negative_affect_1 ~ negative_affect_0 + exposure_mean + (1|team_ID), data = node)
model_slp = lmer(negative_affect_1 ~ negative_affect_0 + exposure_mean + (exposure_mean|team_ID), data = node)

# R Code: view results
summary(model_OLS)
summary(model_int)
summary(model_slp)

