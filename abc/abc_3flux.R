#=================================================
# ABC for the Anderson Model - 3 outflows version 
#=================================================

data.folder = "C:/Users/rfuchs/Documents/These/Oceano/carbon_pump_abc/abc" 
res.folder = file.path(data.folder,"simus/simus_3flows_20201203/") 
code.folder = "C:/Users/rfuchs/Documents/GitHub/CarbonPump/abc"
source(file.path(code.folder, 'model.R'))

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
    d = euc.dist(matrix(pred,1), as.numeric(matrix(stats,1)))
    distances = append(distances, d)
  }
  
  for (i in 1:n){
    idx = nth(distances, i, descending = F, index.return = T)
    best_params[i,] = c(params[idx,])
  }
  best_params
}

stats.draws <- function(params, nb.outflows){
  # Draw samples from the posterior of the summary statistics
  nb_boot = dim(params)[1]
  # ! The number of summary statistics is hardcoded here
  draws = matrix(NA, nb_boot, nb.outflows)
  
  for (i in 1:nb_boot){
    pred = anderson(params[i,])
    draws[i,] = pred
  }
  draws
}

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]


anderson_labels = c('POC', 'DOC', '$\\psi$', '$\\omega_a$', '$\\omega_{FL}$', '$\\alpha$', '$CF_a$', 
                    '$CF_{FL}$')
anderson_labels_txt = c('POC', 'DOC', 'psi', 'wA', 'wFL', 'alpha', 'CF_A', 
                    'CF_FL')

stats_labels = c('Production_NonSinking', 'Production_Sinking', 'Respiration_zoo')

#=========================================
# Launching the model
#=========================================

# Priors and summary statistics
nb_params = length(anderson_labels)
nb_fixed_entries = length(c('poc', 'doc'))
priors <- vector("list", nb_params)

# The first 3 flows are the data of POC and DOC and are ruled out later
priors[[3]] = c("unif",0,1)
priors[[4]] = c("unif",0,0.2)
priors[[5]] = c("unif",0,0.2)
priors[[6]] = c("unif",0,1)
priors[[7]] = c("unif",0,10)
priors[[8]] = c("unif",0,10)

out.flows = read.csv(file.path(data.folder, 'out_flows.csv'))  
names(out.flows) = c('Cruise', 'Station', 'Depths',
                     'POC net', 'DOC net', 'prod non sinking',
                     'prod sinking', 'respi zoo')
sum_stat_obs = out.flows[,c('prod non sinking',
                            'prod sinking', 'respi zoo')]

nb_station = dim(sum_stat_obs)[1]
nb_stats = length(sum_stat_obs)

# Hyper parameters
n = 10000 # 25 minutes for 20000
p = 0.0005

# All station simulations are stored in ABC_rej
ABC_rej = vector("list", nb_station)

for (loc.idx in 1:nb_station){
  
  poc = out.flows[loc.idx,'POC net']
  doc = out.flows[loc.idx,'DOC net']

  cruise = out.flows[loc.idx,'Cruise']
  station = out.flows[loc.idx,'Station']
  station = gsub('/', '-', station)
  print(paste(c(station,cruise)))
  
  # The tree first flux are the data extracted for each site
  priors[[1]] = c("unif", poc, poc)
  priors[[2]] = c("unif", doc, doc)

  start = Sys.time()
  # Est-ce que les stats sont bien normalisées ? :)
  # A checker: pk renvoie 6 paramètres et pas 9 ?
  ABC_rej[[loc.idx]] <- ABC_rejection(model = anderson, prior = priors, 
                                      nb_simul=n, 
                                      summary_stat_target = c(sum_stat_obs[loc.idx,]), 
                                      tol=p)
  end  = Sys.time()
  print(end - start)

  
  #=====================================
  # Plot the stats posterior densities
  #=====================================
  
  # As the first 3 params are constant there are 
  # not returned by ABC, so we re-add them.
  params = ABC_rej[[loc.idx]]$param
  poc.doc = t(replicate(n * p, c(poc, doc)))
  params = cbind(poc.doc, params)
  
  # Compute the three outflows predicted by the model for 
  # the p "best" simulations
  stats.posterior = stats.draws(params, nb.outflows = nb_stats)

  #===============================
  # Joint density
  #===============================
  
  # Determine the best joint combinaison of parameters
  nb_best_points = 1
  best_joint_params = best.comb(params, c(sum_stat_obs[loc.idx,]), nb_best_points)
  best_pred = anderson(best_joint_params[1,])
  
  # Plot outflows and save the figure
  png(file=file.path(res.folder, 'outflows',  paste(cruise, station, ".png")),
      width=800, height=550)
  par(mfrow=c(2,2))
  
  for (k in 1:nb_stats){
    
    max.x.displayed = max(sum_stat_obs[loc.idx, k] + 5, max(stats.posterior[,k]))
    hist(stats.posterior[,k], main = stats_labels[k], 
         breaks = 20, xlab = "Predicted Value",
         xlim = c(0, max.x.displayed))
    
    abline(v = sum_stat_obs[loc.idx, k], col = 'green', main = "") 
    abline(v = best_pred[k], col = 'blue', main = "") 
    legend(1, 1, legend=c("True flow", "Best pred"),
           col=c("green", "blue"), lty=1:2, cex=0.65)
  }
  mtext(paste(cruise, station), side = 3, line = -1.5, outer = TRUE)
  dev.off()
  
  
  
  # Plotting the marginal densities along with the best params
  png(file=file.path(res.folder, 'params_plots',  paste(cruise, station, ".png")),
      width=800, height=550)
  par(mfrow=c(2,3)) 
  
  for (k in (nb_fixed_entries + 1):nb_params){
    hist(params[,k], main = TeX(anderson_labels[k]),
         xlab = 'Predicted value',
         breaks = 20)
    abline(v = best_joint_params[1,k], col = 'blue')
    legend(0.1, 0.5, legend=c("Best pred"),
           col=c("blue"), lty=1:2, cex=0.65)
    }
  mtext(paste(cruise, station), side = 3, line = -1.5, outer = TRUE)
  dev.off()
  
  # save the params drawn
  n_best = matrix(params, n * p, nb_params)
  write.table(n_best[, (nb_fixed_entries + 1):nb_params], 
            col.names = anderson_labels_txt[(nb_fixed_entries + 1):nb_params],
            file.path(res.folder, 'params', 'n_best',
            paste(cruise, station, ".csv")),
            row.names = F)
  
  one.best = matrix(best_joint_params[,(nb_fixed_entries + 1):nb_params], 
                    nb_best_points, nb_params - nb_fixed_entries )
  
  write.table(one.best, 
              col.names = anderson_labels_txt[(nb_fixed_entries + 1):nb_params],
              file.path(res.folder, 'params', '1_best', 
              paste(cruise, station, ".csv")),
              row.names = F)
}

# Read all best params files and create a summary csv file:
files = list.files(file.path(res.folder, 'params', '1_best'))
files = files[files != 'all.csv']

best.params.all = as.data.frame(matrix(NA, length(files), nb_params - nb_fixed_entries))
for (f.idx in 1:length(files)) {
  current.params = read.table(file.path(res.folder, 'params', '1_best', files[f.idx]),
                       header = T)
  best.params.all[f.idx,] = c(current.params[1,])
}

filenames = gsub(' .csv', '', files)
write.table(best.params.all, 
            col.names = anderson_labels_txt[(nb_fixed_entries + 1):nb_params],
            row.names = filenames,
            file.path(res.folder, 'params', '1_best', 
            "all.csv"), sep = ',')


# Test d'affichage d'histogramme 2D
library(ggplot2)

cfs = data.frame(n_best[,7:8])
names(cfs) = c('x','y')

ggplot(cfs, aes(x=x, y=y) ) +
  geom_bin2d() +
  theme_bw()
