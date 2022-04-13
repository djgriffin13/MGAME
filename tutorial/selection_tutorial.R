# Install and load library
devtools::install_github("djgriffin13/isNet")
library(isNet)

# Lode data the path needs to go to the data
data = readr::read_csv("./tutorial/data/selection_data.csv")

# Calculate distance
data$na_dist = sqrt((data$na_source - data$na_target) ^ 2)

# Run model
fit = ameMG(
  data,
  Y = "communication_count_t1",
  Xego = c("na_source"),
  Xalter = "na_target",
  Xdyad = c("communication_count_t0", "na_dist"),
  group = "team_ID",
  R = 2,
  nscan = 300,
  family = "nrm"
)

# Get results
summary(fit)
