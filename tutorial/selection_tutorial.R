# Install and Load Library
# If devtools is not already installed
# it allows you to install MGAME from github directly
install.packages("devtools")
devtools::install_github("djgriffin13/MGAME", upgrade = "never")
library(MGAME)

# Load Data
# Must set path to local directory where the data is stored
# Click: Session > Set Working Directory > Choose Directory...
data = readr::read_csv("./data/selection_data.csv")

# Calculate Relational Covariates
data$negative_affect_sim_t0 = 5 - abs(data$negative_affect_ego_t0 - data$negative_affect_alter_t0)

# Run Model
fit = mgame(
   data,
   Y = "info_share_t1",
   Xego = c("negative_affect_ego_t0"),
   Xalter = c("negative_affect_alter_t0"),
   Xdyad = c("info_share_t0", "negative_affect_sim_t0"),
   group = "team_ID",
   R = 2,
   nscan = 10000,
   family = "nrm",
   group_standard = c("negative_affect_ego_t0","negative_affect_alter_t0","negative_affect_sim_t0"),
   grand_standard = c("info_share_t0","info_share_t1"),
)

# Get Results
summary(fit)
