#====================================================
# Optimisation : Anderson Model - 3 outflows version 
#====================================================

#install.packages('dfoptim')
library(dfoptim)

data.folder = "C:/Users/rfuchs/Documents/These/Oceano/carbon_pump_abc/opt/" 
res.folder = "C:/Users/rfuchs/Documents/These/Oceano/carbon_pump_abc/opt/simus/20211401/"
code.folder = "C:/Users/rfuchs/Documents/GitHub/CarbonPump/optimisation"
source(file.path(code.folder, 'model.R'))

#======================================
# Utilities definition
#======================================

euc.dist <- function(x1, x2, w = 1) sqrt(rowSums((w *(x1 - x2)) ^ 2))

anderson_labels_txt = c('POC', 'DOC', 'psi', 'wA', 'wFL', 'alpha', 'CF_A', 
                        'CF_FL')

pred.error <- function(X){
  models.params = X[1:8]
  out.flows = as.numeric(X[9:12])
  pred = anderson(as.numeric(models.params))  
  w = 1 / as.numeric(out.flows)
  d = euc.dist(matrix(pred,1), matrix(out.flows,1), w)
  d
}

#=========================================
# Launching the model
#=========================================

# Priors and summary statistics
nb_params = length(anderson_labels_txt)
nb_fixed_entries = length(c('poc', 'doc'))

flows = read.csv(file.path(data.folder, 'out_flows_4flows_PAP.csv'), sep = ';')  
names(flows) = c('Cruise', 'Station', 'Depths',
                     'POC net', 'DOC net', 'prod non sinking',
                     'prod sinking', 'resp sinking', 'resp zoo')

out.flows = flows[,c('prod non sinking', 'prod sinking', 
                     'resp sinking', 'resp zoo')]

nb_stations = dim(out.flows)[1]
nb_stats = length(out.flows)

# All station simulations are stored in ABC_rej
labels = anderson_labels_txt[(nb_fixed_entries + 1):nb_params]

res = data.frame(matrix(NA, nb_stations, 
                nb_params - nb_fixed_entries + nb_stats))

colnames(res) = c(labels, c('prod non sinking','prod sinking',
                            'resp sinking', 'resp zoo'))

for (loc.idx in 1:(nb_station + 1)){
  epsilon = 1E-8
  
  poc = flows[loc.idx,'POC net']
  doc = flows[loc.idx,'DOC net']
  
  cruise = flows[loc.idx,'Cruise']
  station = flows[loc.idx,'Station']
  station = gsub('/', '-', station)
  print(paste(c(station,cruise)))
  
  # The initial parameters values to launch the algorithm
  # Check starting values of the PGEs in Giering/literature
  starting.values = c(c(poc, doc, 0.75, 0.02, 0.08, 0.5, 0.44, 0.44), out.flows[loc.idx,]) 
  starting.values = as.numeric(starting.values)
  # Lower bounds contraints for each parameter
  lb = c(c(poc - epsilon, doc - epsilon), rep(0, nb_params - nb_fixed_entries), as.numeric(out.flows[loc.idx,]) - epsilon) 
  lb = as.numeric(lb)
  ub = c(c(poc + epsilon, doc + epsilon, 1, 0.1, 0.1, 1, 6, 6), as.numeric(out.flows[loc.idx,]) + epsilon) # Lower bounds contraints for each parameter
  ub = as.numeric(ub)
  
  start = Sys.time()
  res.station = nmkb(par = starting.values, fn = pred.error, lower=lb, upper=ub)
  end  = Sys.time()
  print(end - start)
  
  # Compute the error existing between flows
  preds = anderson(res.station$par[1:8])
  flows.error = as.numeric(out.flows[loc.idx,]) - preds
  print(paste('flow error', flows.error))

  # save the params drawn
  best.params.values = res.station$par[3:8]
  best.params.values = matrix(best.params.values, 1, nb_params - nb_fixed_entries)

  res[loc.idx,] = as.numeric(c(best.params.values, c(flows.error)))
    
  write.table(best.params.values, 
              col.names = T,
              file.path(res.folder, 'params', 
              paste(cruise, station, ".csv")),
              row.names = F)
}

files = list.files(file.path(res.folder, 'params'))
files = files[files != 'all.csv']
filenames = gsub(' .csv', '', files)

write.table(res, col.names = T,
            file.path(res.folder,
            'params', "all.csv"),
            row.names = filenames)
print(res)

  