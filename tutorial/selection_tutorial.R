# Run install.packages("devtools")
# if it is not already installed
# this allows you to install from github directly
install.packages("devtools")

# Install and Load MGAME package from github
devtools::install_github("djgriffin13/MGAME")
library(MGAME)

# Load Data
# Mut set path to local directory where the data is stored
data = readr::read_csv("./tutorial/data/selection_data.csv")

# Calculate Dissimilarity
data$na_sim = 5 - abs(data$na_ego - data$na_alter)

# Run Model
fit = mgame(
  data,
  Y = "info_share_t1",
  Xego = c("na_ego"),
  Xalter = "na_alter",
  Xdyad = c("info_share_t0", "na_sim"),
  group = "team_ID",
  R = 2,
  nscan = 300,
  family = "nrm"
)

# Get Rresults
summary(fit)
