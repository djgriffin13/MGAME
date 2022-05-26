# Import Packages
install.packages("tidyverse")
library(tidyverse)

# Import Data
edge = read_csv("./data/edge_influence.csv")
node = read_csv("./data/node_influence.csv")

# Merge Data Sets to the Dyadic Level
merged_data = left_join(edge,node, by = c("alter_ID" = "individual_ID"))

# Calculate the Dyadic Exposure
merged_data = merged_data %>% mutate(exposure_dyad = information_receipt_t0 * negative_affect_t0)

# Calculate the Exposure Term
total_exposure_table = merged_data %>% group_by(ego_ID) %>% summarize(exposure_t0 = sum(exposure_dyad,na.rm = TRUE))

# Merge Data Sets to the Individual Level
final_data = left_join(node,total_exposure_table, by = c("individual_ID" = "ego_ID"))

# Influence Model in OLS Regression (both unstandardized and standardized)
model_OLS = lm(negative_affect_t1 ~ negative_affect_t0 + exposure_t0, data = final_data)
model_OLS_std = lm(scale(negative_affect_t1) ~ scale(negative_affect_t0) + scale(exposure_t0),
                   data = final_data)

# View Results
summary(model_OLS) # Get unstandardized regression Coefficents
summary(model_OLS_std) # Get Standardized Regression Coefficents
confint.default(model_OLS_std) #Get confidence Intervals
coef(model_OLS_std)["exposure_t0_std"]/coef(model_OLS_std)["negative_affect_t0_std"]

