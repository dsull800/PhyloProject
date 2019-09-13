require("RPANDA")
require("diversitree")
require("castor")
require("phybase")
require("ggtree")
require("treeAGG")
require("ape")

# define an HBD model with exponentially decreasing speciation/extinction rates
Ntips = 1000
beta = 0.01 # exponential decay rate of lambda over time
oldest_age= 10
age_grid = seq(from=0,to=oldest_age,by=0.1) # choose a sufficiently fine age grid
lambda = 1*exp(beta*age_grid) # define lambda on the age grid
mu = 0.2*lambda # assume similarly shaped but smaller mu
# simulate deterministic HBD model
simulation = simulate_deterministic_hbd(Ntips,
                                        oldest_age,rho = 0.5,
                                        age_grid = age_grid,
                                        lambda = lambda,
                                        mu = mu)
# plot deterministic LTT
plot( x = simulation$ages, y = simulation$LTT, type='l',
      main='dLTT', xlab='age', ylab='lineages')