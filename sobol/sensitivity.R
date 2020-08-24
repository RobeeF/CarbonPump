#############
#* Sensitivity
#############

setwd("C:/Users/rfuchs/Documents/GitHub/CarbonPump/sobol")
source('model.R')
source('model_seq.R')

#install.packages("TSP", dependencies = TRUE)
#install.packages('RANN')
#install.packages('sensitivity')

library(sensitivity)
library(ggplot2)
library(latex2exp)


anderson_labels = c("psi", 'omega_a', 'omega_FL', 'alpha', 'phi_v', 
                    'beta_v', 'K_v', 'phi_v', 'beta_v', 'K_v',
                    'phi_z', 'beta_z', 'lambda_z', 'K_z', 'phi_h',
                    'beta_H', 'lambda_h', 'K_h', 'zeta', 'zi2')
TeX(anderson_labels[1])


#==================================================
# Sobol knn
#==================================================
N <- 50000
k <- 20
X <- data.frame(matrix(runif(k * N), nrow = N))

stats_labels = c('Production_NonSinking', 'Production_Sinking', 'Respiration_Sinking', 'Respiration_zoo')


order.one.indices = data.frame(matrix(0, 4, 20))
colnames(order.one.indices) <- anderson_labels

total.indices = matrix(0, 4, 20)
colnames(total.indices) <- anderson_labels


for (flow in 1:4){
  cat('flow number', flow, '\n')
  start = Sys.time()
  x <- sobolshap_knn(model = anderson, X = X, U = 1,
                     method = "rank", flux_out_nb = flow)

  
  colnames(x$S) <- anderson_labels
  print(ggplot(x)  + ggtitle(paste('1st order Sobol indices:', stats_labels[flow])))
  order.one.indices[flow,] = x$S
  
  x2 <- sobolshap_knn(model = NULL, X = X, U = 0, method = "knn", n.knn = 5)
  tot.indices = tell(x2,x$y)
  
  colnames(x2$S) <- anderson_labels
  total.indices[flow,] = x2$S
  print(ggplot(x2)  + ggtitle(paste('Total Sobol indices:', stats_labels[flow])))
  
  end = Sys.time()
  print(end-start)
  
}

write.csv(order.one.indices, '1st_order_Sobol.csv')
write.csv(total.indices, 'total_Sobol.csv')


for (flow in 1:4){
  x$S = order.one.indices[flow,]
  print(ggplot(x)  + ggtitle(paste('First order Sobol indices:', stats_labels[flow])))
}

for (flow in 1:4){
  plot(total.indices[flow,], main = paste('Total Sobol indices:', stats_labels[flow]),
       ylab = 'Total Sobol indices')
  axis(1, at=1:20, labels = anderson_labels)
}


#===========================================
# Legacy code
#===========================================

# Start with the first order indices
start = Sys.time()
x <- sobolshap_knn(model = anderson, X = X, U = 1,
                   method = "rank", flux_out_nb = 1)
end = Sys.time()
print(end-start)
# 2000 particles = 9s, 5000pts = 22s, 10000pts = 1,6mins, 20 000pts = 4,6mins, 40 000 pts = 43mins
# 50 000 pts 66 minutes

colnames(x$S) <- anderson_labels
ggplot(x)
xx <- as.data.frame(t(as.data.frame(x$S)))

plot <- ggplot(xx, aes(anderson_labels))


# Then the total indices
x2 <- sobolshap_knn(model = NULL, X = X, U = 0, method = "knn", n.knn = 5)
tot.indices = tell(x2,x$y)
ggplot(x2)  + scale_x_discrete(labels = c('X1' = anderson_labels[1]))



