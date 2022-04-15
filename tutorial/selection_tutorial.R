# Install and Load Library
# Note only run install.packages("devtools")
# if it is not already installed
install.packages("devtools")
devtools::install_github("djgriffin13/isNet")
library(isNet)

# Load Data
# Set path to local directory where the data is stored
data = readr::read_csv("./tutorial/data/selection_data.csv")

# Calculate Dissimilarity
data$na_dissim = sqrt((data$na_source - data$na_target) ^ 2)

# Run Model
fit = ameMG(
  data,
  Y = "communication_count_t1",
  Xego = c("na_source"),
  Xalter = "na_target",
  Xdyad = c("communication_count_t0", "na_dissim"),
  group = "team_ID",
  R = 2,
  nscan = 300,
  family = "nrm"
)

# Get Rresults
summary(fit)
