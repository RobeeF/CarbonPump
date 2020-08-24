setwd("C:/Users/rfuchs/Documents/GitHub/CarbonPump/sobol")
source('model_seq.R')
source('model.R')

#install.packages('sensobol')
#install.packages('hexbin')

library(sensobol)
library(hexbin)
library(ggplot2) # To plot the Sobol' bootstrap replicas

#===========================
# Check results concordance
#===========================

# To check before launching the following
N <- 2
k = 20
params <- paste("X", 1:k, sep = "")

## Create sample matrix to compute first and total-order indices:
mat <- sobol_matrices(N, k)
length(anderson(mat))

res = anderson(mat)

for (j in 1:44){
  res_old = anderson_seq(mat[j,])

  print(res[j,] - res_old)
}

# No error for the parallelized version

#===========================
# Benchmark running time:
#===========================

para = rep(0, 5)
non_para = rep(0, 5)

N = c(40, 100, 200, 400, 1000)

for (i in 1:5){
  cat('Number of points:', N[i], '\n')
  mat <- sobol_matrices(N[i], k)
  n_comb = dim(mat)[1]
  
  start = Sys.time()
  res = anderson(mat)
  end = Sys.time()
  para[i] = ((end - start))
  
  
  start = Sys.time()
  for (j in 1:n_comb){
    res = anderson_seq(mat[j,])
  }
  end = Sys.time()
  non_para[i] = ((end - start))
}

plot(para, type ='l', col = 'green')
lines(non_para , type ='l', col ='red')

# Parallelized is slower for high N...

#================================
# Compute Sobol indices
#================================

# All the original code is available at :
# https://cran.r-project.org/web/packages/sensobol/vignettes/sensobol.html

# Define the parameters

N <- 5000 # Sample size
k <- 20 # Number of parameters
params <- paste("X", 1:k, sep = "") # Vector with the name of the model inputs
R <- 100 # Number of bootstrap replicas

# Create the Sobol' matrices
A <- sobol_matrices(n = N, 
                    k = k) 

# Create sample matrix to compute first and total-order indices:
n_comb = dim(A)[1]

'''
start = Sys.time()
for (j in 1:n_comb){
  Y[j,] = anderson_seq(A[j,])
}
end = Sys.time()
print(end-start)
'''

start = Sys.time()
Y = anderson(A)
end = Sys.time()
print(end-start)

# Perform the analysis for each variable separately
Yk = Y[,1] # Take 1 for instance

#Y <- anderson(A)

# Plot the repartition of Yk for all bootstrap draws
plot_uncertainty(Yk, n = N)

# Draw Yk = f(inputs_i) for all i in [1,20]
plot_scatter(x = A, 
             Y = Yk, 
             n = N, 
             params = params)

# Compute the Sobol indices of first order and total (can compute for order 2 and 3 if needed)
dt <- sobol_indices(Y = Yk, 
                    params = params, 
                    type = "saltelli", 
                    R = R,
                    n = N, 
                    second = FALSE, 
                    third = FALSE)

# Extract the bootstrap samples
b.rep <- sobol_replicas(dt = dt, k = k)


# Plot of the Sobol indices (order 1 and total)
ggplot2::ggplot(b.rep, aes(value)) +
  geom_histogram() +
  labs(x = "Y",
       y = "Count") +
  facet_wrap(parameters~variable, 
             scales = "free") +
  labs(x = "Variance", 
       y = "Count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent",
                                         color = NA),
        legend.key = element_rect(fill = "transparent",
                                  color = NA)) 

# Compute confidence intervals 
# Take gaussian interval
dt.ci <- sobol_ci(dt, 
                  params = params, 
                  type = "norm", 
                  conf = 0.95, 
                  second = FALSE, 
                  third = FALSE) 


dt.dummy <- sobol_dummy(Y = Yk, 
                        params = params, 
                        R = R, 
                        n = N)

dt.dummy.ci <- sobol_ci_dummy(dt.dummy, 
                              type = "norm", 
                              conf = 0.95)

plot_sobol(dt.ci, dummy = dt.dummy.ci, type = 1)


#===============================================
# Automatic detection of influencial variables
#===============================================

signif_coef1 = list()
signif_coef2 = list()

var_thr= 0.1

for (j in 1:85){
  cat(j,'th flow \n')
  Yk = Y[,j]
  
  dt <- sobol_indices(Y = Yk, 
                      params = params, 
                      type = "saltelli", 
                      R = R,
                      n = N, 
                      second = FALSE, 
                      third = FALSE)
  
  b.rep <- sobol_replicas(dt = dt, k = k)
  
  dt.ci <- sobol_ci(dt, 
                    params = params, 
                    type = "norm", 
                    conf = 0.95, 
                    second = FALSE, 
                    third = FALSE) 
  
  
  dt.dummy <- sobol_dummy(Y = Yk, 
                          params = params, 
                          R = R, 
                          n = N)
  
  dt.dummy.ci <- sobol_ci_dummy(dt.dummy, 
                                type = "norm", 
                                conf = 0.95)

  # We select only over 1st order coefficients
  order1 = dt.ci[dt.ci$sensitivity == 'Si']
  order1.dummy = dt.dummy.ci[dt.dummy.ci$sensitivity == 'Si']
  
  # Selecting using confidence intervals
  signif_idx = (order1$original <=  order1.dummy$low.ci) | (order1$original >=  order1.dummy$high.ci)
  signif_coef1 = append(signif_coef1, c(order1[signif_idx,'parameters']))
  
  # Selecting using using threshold
  signif_idx = order1$original >= var_thr
  
  signif_coef2 = append(signif_coef2, c(order1[signif_idx,'parameters']))
  
  print('Total indices sum:')
  print(sum(dt.ci[dt.ci$sensitivity == 'STi']$original))
  print('1st order indices sum:')
  print(sum(dt.ci[dt.ci$sensitivity == 'Si']$original))
  
}
