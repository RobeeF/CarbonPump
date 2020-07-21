setwd("C:/Users/rfuchs/Documents/These/Oceano/carbon_pump_abc/abc")
source('model.R')

#install.packages("EasyABC") 
#install.packages("mnormt")
#install.packages("Rfast")
#install.packages('latex2exp')

library(latex2exp)
library(EasyABC)
library(Rfast)

#======================================
# Utilities definition
#======================================

dmode <- function(x, precision = 512) {
  den <- density(x, kernel = c("gaussian"), n = precision)
  ( den$x[den$y == max(den$y)] )   
} 


euc.dist <- function(x1, x2) sqrt(rowSums((x1 - x2) ^ 2))


best.comb <- function(params, stats, n){
  # Return the n closest combinaisons of parameters
  
  nb_boot = dim(params)[1]
  nb_params = dim(params)[2]
  distances = c()
  
  best_params = matrix(NA, n, nb_params)
  
  for (i in 1:nb_boot){
    pred = anderson(params[i,])
    d = euc.dist(matrix(pred,1), matrix(stats,1))
    distances = append(distances, d)
  }
  
  for (i in 1:n){
    idx = nth(distances, i, descending = F, index.return = T)
    best_params[i,] = c(params[idx,])
  }
  best_params
}

stats.draws <- function(params){
  # Draw samples from the posterior of the summary statistics
  nb_boot = dim(params)[1]
  # ! The number of summary statistics is hardcoded here
  draws = matrix(NA, nb_boot, 4)
  
  for (i in 1:nb_boot){
    pred = anderson(params[i,])
    draws[i,] = pred
  }
  draws
}

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]


anderson_labels = c("$\\psi$", '$\\omega_a$', '$\\omega_{FL}$', '$\\alpha$', '$\\phi_v$', 
                    '$\\beta_v$', '$K_v$', '$\\phi_v$', '$\\beta_v$', '$K_v$',
                    '$\\phi_z$', '$\\beta_z$', '$\\lambda_z$', '$K_z$', '$\\phi_h$',
                    '$\\beta_H$', '$\\lambda_h$', '$K_h$', '$\\zeta$', 'zi2')

stats_labels = c('Respiration attached', 'Respiration zoo', 'Production free living', 'Production attached')


#=========================================
# Launching the model
#=========================================

# Priors and summary statistics
n_params = 20
priors <- vector("list", n_params)
for (i in 1:n_params) priors[[i]] <- c("unif",0,1)

sum_stat_obs = c(60.19, 1.62, 17.75, 40.46)
#sum_stat_obs = c(17.75345, 40.46, 58.067, 0.28)
nb_stats = length(sum_stat_obs)

# Hyper parameters
n = 50000 # 25 minutes for 20000
p = 500/n

# Running the algorithm
start = Sys.time()
ABC_rej <- ABC_rejection(model = anderson, prior = priors, nb_simul=n, summary_stat_target = sum_stat_obs, tol=p)
end  = Sys.time()
print(end - start)

#=====================================
# Plot the stats posterior densities
#=====================================
params = ABC_rej$param

stats.posterior = stats.draws(params)
par(mfrow=c(2,2))

for (k in 1:nb_stats){
  plot(density(stats.posterior[,k]), main = stats_labels[k])
  abline(v = sum_stat_obs[k], col = 'blue')  
  
}


#===============================
# Marginal densities (useless)
#===============================

par(mfrow=c(4,5)) 

# Plotting the marginal densities
for (k in 1:n_params){
plot(density(ABC_rej$param[,k]), main = paste0('X', k))
}

# Storing the more likely parameters according to marginal densities
best_marg_params = rep(0, n_params)
for (k in 1:n_params) best_marg_params[k] <- dmode(ABC_rej$param[,k], precision = 1e4)

# Model output for this set of params
print('True respirations/productions:')
print(sum_stat_obs)
print('Predicted respirations/productions:')
pred = anderson(best_marg_params)
print(pred)
print('Distance:')
print(euc.dist(matrix(pred,1), matrix(sum_stat_obs,1)))



#===============================
# Joint density
#===============================

# Determine the best joint combinaison of parameters
nb_best_points = 5
best_joint_params = best.comb(ABC_rej$param, sum_stat_obs, nb_best_points)

# Model output for this set of params
print('True respirations/productions:')
print(sum_stat_obs)
print('Predicted respirations/productions:')
pred = anderson(best_joint_params[1,])
print(pred)
print('Distance:')
print(euc.dist(matrix(pred,1), matrix(sum_stat_obs,1)))


# Choose colors
cols = sample(color, nb_best_points)

# Plotting the marginal densities along with the best params
par(mfrow=c(4,5)) 


for (k in 1:n_params){
  plot(0, 0, type = 'n', main = TeX(anderson_labels[k]))
  
  #plot(density(ABC_rej$param[,k]), main = TeX(anderson_labels[k]))

  for (bp in 1:3){
  abline(v = best_joint_params[bp], col = cols[bp])  
  }
}


# save the params drawn
write.csv(matrix(ABC_rej$param, 500, 20), 'params.csv')
write.csv(matrix(best_joint_params[1,], 1, 20), 'best_params.csv')

