# Install needed packages
install.packages("tidyverse")

# R Code: Load data
library(tidyverse)
edge = read_csv("./tutorial/data/edge_influence.csv")
node = read_csv("./tutorial/data/node_influence.csv")


# R Code: Join node data to edge by source id
full_edge_data = left_join(edge,node, by = c("ego_ID" = "node_ID"))

# R Code: Calculate Dyadic exposure
full_edge_data = full_edge_data %>% mutate(exposure = information_receipt * negative_affect_t0)

# Calculate mean exposure 
total_exposure_tbl = full_edge_data %>% group_by(alter_ID) %>% summarize(total_exposure = sum(exposure,na.rm = TRUE))

# Include mean exposure in the node's data and handle cases where there was no exposure
node = left_join(node,total_exposure_tbl, by = c("node_ID" = "alter_ID"))

# R Code: Run OLS analysis
model_OLS = lm(negative_affect_t1 ~ negative_affect_t0 + total_exposure, data = node)

# R Code: view results
summary(model_OLS)


