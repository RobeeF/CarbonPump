#############
#* Sensitivity
#############

setwd("C:/Users/rfuchs/Documents/These/Oceano/carbon_pump_abc/sobol")
source('model.R')
source('old_model.R')

#install.packages("TSP", dependencies = TRUE)
#install.packages('RANN')
#install.packages('sensitivity')

library(sensitivity)
library(ggplot2)

#==================================================
# Sobol knn
#==================================================
N <- 5000 
k <- 20
X <- data.frame(matrix(runif(k * N), nrow = N))

# Start with the first order indices
start = Sys.time()
x <- sobolshap_knn(model = anderson, X = X, U = 1,
                   method = "rank", flux_out_nb = 1)
end = Sys.time()
print(end-start)
# 2000 particles = 9s, 5000pts = 22s, 10000pts = 100s, 20 000pts = 281s

order1.indices <- x$S
ggplot(x)


# Then the total indices
x2 <- sobolshap_knn(model = NULL, X = X, U = 0, method = "knn", n.knn = 5)
tot.indices = tell(x2,x$y)
ggplot(x2)
