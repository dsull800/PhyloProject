require("RPANDA")
require("diversitree")
require("castor")
require("phybase")
require("ggtree")
require("treeAGG")
require("ape")

# t <- 0:10  # time
# real_temp_data=data.frame(t,t)
# f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
# f.mu<-function(t,x,y){0}
# lamb_par<-c(0.2, 0.01)
# mu_par<-c()
# result_exp <- sim_env_bd(real_temp_data,f.lamb,f.mu,lamb_par,mu_par,time.stop=100)
# plot.phylo(result_exp[[1]])

#actually instead of using code above use castor code generate_random_tree

## Not run:
# Generate a random tree with exponentially varying lambda & mu
Ntips = 20
rho = 0.5 # sampling fraction
time_grid = seq(from=0, to=100, by=0.01)
lambdas = 2*exp(0.1*time_grid)
mus = 0*time_grid
sim = generate_random_tree( parameters = list(rarefaction=rho),
                            max_tips = Ntips/rho,
                            coalescent = FALSE,
                            added_rates_times = time_grid,
                            added_birth_rates_pc = lambdas,added_death_rates_pc = mus)
bigtree = sim$tree
root_age = castor::get_tree_span(bigtree)$max_distance
cat(sprintf("Tree has %d tips, spans %g Myr\n",length(bigtree$tip.label),root_age))
bigtree$tip.label=1:length(bigtree$tip.label)
plot.phylo(bigtree)

#this code below gets the rootnode of the phylogeny but not actually
#mat=matTree(tree=result_exp[[1]
 
# sim.coaltree.sp(rootnode=mat[length(mat)],nodematrix=mat,nspecies=result_exp[["nblineages"]][length(result_exp[["nblineages"]])],seq=rep(1,result_exp[["nblineages"]][length(result_exp[["nblineages"]])]))

newickformat=write_tree(bigtree)

# write.tree(phy=result_exp[[1]],file="newick1.txt")
# realtree=read.tree(file = "newick1")

#Be careful, the newickformat string ends in \n strsplit(newickformat,"\n")
# fileName <- 'newick1.txt'
# newickformat=readChar(fileName, file.info(fileName)$size)

# newickformat=read.tree.string("newick1.tre",format="phylip")

nodematrix=read.tree.nodes(newickformat)

#I can set the population sizes by changing the entries of the 5th value of the matrix, 
#use a gamma distributed value or something, maybe.
for(i in 1:nrow(nodematrix$nodes)){
  nodematrix$nodes[i,5]=rgamma(n=1,shape=2,rate=1)
}

# nspecies=tail(result_exp[[3]],n=1)

nspecies=length(bigtree$tip.label)

specnames=1:20
  #as.integer(bigtree$tip.label)

#nodematrix$nodes[nrow(nodematrix$nodes)-1,1]
rootnode=nrow(nodematrix$nodes)

genetreestuff=sim.coaltree.sp(rootnode,nodematrix$nodes,nspecies,seq=rep(1,nspecies),name=specnames)

genetreeheight=genetreestuff$height

realgenetree=genetreestuff$node

genetreegt=genetreestuff$gt

# plot.phylo(genetreegt)

thetree=read.tree(text=genetreegt)

plot.phylo(thetree)

# Define a parametric HBD congruence class, with exponentially varying PDR
# The model thus has 3 parameters
PDR_function = function(ages,params){
  return(params['A']*exp(-params['B']*ages));
}
rholambda0_function = function(params){
  return(params['rholambda0'])
}
# Define an age grid on which lambda_function & mu_function shall be evaluated
# Should be sufficiently fine to capture the variation in the PDR
age_grid = seq(from=0,to=100,by=0.01)
# Perform fitting
# Lets suppose extinction rates are already known
cat(sprintf("Fitting class to tree..\n"))
fit = fit_hbd_pdr_parametric( thetree,
                              param_values = c(A=NA, B=NA, rholambda0=NA),
                              param_guess = c(1,0,1),
                              param_min = c(-10,-10,0),
                              param_max = c(10,10,10),
                              param_scale = 1, # all params are in the order of 1
                              PDR = PDR_function,
                              rholambda0 = rholambda0_function,
                              age_grid = age_grid,
                              Ntrials = 10, # perform 10 fitting trials
                              Nthreads = 2, # use 2 CPUs
                              max_model_runtime = 1, # limit model evaluation to 1 second
                              fit_control = list(rel.tol=1e-6))
if(!fit$success){
  cat(sprintf("ERROR: Fitting failed: %s\n",fit$error))
}else{
  cat(sprintf("Fitting succeeded:\nLoglikelihood=%g\n",fit$loglikelihood))
  print(fit)
}