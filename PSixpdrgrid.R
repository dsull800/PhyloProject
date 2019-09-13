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
sim = castor::generate_random_tree( parameters = list(rarefaction=rho),
                            max_tips = Ntips/rho,
                            coalescent = FALSE,
                            added_rates_times = time_grid,
                            added_birth_rates_pc = lambdas,added_death_rates_pc = mus)
bigtree = sim$tree
root_age = castor::get_tree_span(bigtree)$max_distance
cat(sprintf("Tree has %d tips, spans %g Myr\n",length(bigtree$tip.label),root_age))
bigtree$tip.label=1:length(bigtree$tip.label)
ape::plot.phylo(bigtree)
title("speciestree")
#this code below gets the rootnode of the phylogeny but not actually
#mat=matTree(tree=result_exp[[1]

# sim.coaltree.sp(rootnode=mat[length(mat)],nodematrix=mat,nspecies=result_exp[["nblineages"]][length(result_exp[["nblineages"]])],seq=rep(1,result_exp[["nblineages"]][length(result_exp[["nblineages"]])]))

newickformat=castor::write_tree(bigtree)

# write.tree(phy=result_exp[[1]],file="newick1.txt")
# realtree=read.tree(file = "newick1")

#Be careful, the newickformat string ends in \n strsplit(newickformat,"\n")
# fileName <- 'newick1.txt'
# newickformat=readChar(fileName, file.info(fileName)$size)

# newickformat=read.tree.string("newick1.tre",format="phylip")

nodematrix=phybase::read.tree.nodes(newickformat)

#I can set the population sizes by changing the entries of the 5th value of the matrix, 
#use a gamma distributed value or something, maybe.
for(i in 1:nrow(nodematrix$nodes)){
  nodematrix$nodes[i,5]=rgamma(n=1,shape=2,rate=1)
}

# nspecies=tail(result_exp[[3]],n=1)

nspecies=length(bigtree$tip.label)

specnames=1:length(bigtree$tip.label)
#as.integer(bigtree$tip.label)

#nodematrix$nodes[nrow(nodematrix$nodes)-1,1]
rootnode=nrow(nodematrix$nodes)

genetreestuff=phybase::sim.coaltree.sp(rootnode,nodematrix$nodes,nspecies,seq=rep(1,nspecies),
                                       name=specnames)

genetreeheight=genetreestuff$height

realgenetree=genetreestuff$node

genetreegt=genetreestuff$gt

# plot.phylo(genetreegt)

thetree=ape::read.tree(text=genetreegt)

ape::plot.phylo(thetree)
title("genetree")

# calculate true PDR
lambda_slopes = diff(lambdas)/diff(time_grid);
lambda_slopes = c(lambda_slopes[1],lambda_slopes)
PDRs = lambdas - mus - (lambda_slopes/lambdas)
# Fit PDR on grid
Ngrid = 10
age_grid = seq(0,genetreeheight,length.out=Ngrid)
fit = castor::fit_hbd_pdr_on_grid(thetree,
                          age_grid=age_grid,
                          min_PDR = -100,
                          max_PDR = +100,
                          condition = "crown",
                          Ntrials = 10,# perform 10 fitting trials
                          Nthreads = 10,# use two CPUs
                          max_model_runtime = 1) # limit model evaluation to 1 second
if(!fit$success){
  cat(sprintf("ERROR: Fitting failed: %s\n",fit$error))
}else{
  cat(sprintf("Fitting succeeded:\nLoglikelihood=%g\n",fit$loglikelihood))
  # plot fitted & true PDR
  plot( x = fit$age_grid,
        y = fit$fitted_PDR,
        main = 'Fitted & true PDR',
        xlab = 'age',
        ylab = 'PDR',
        type = 'b',
        col = 'red',
        xlim = c(root_age,0))
  lines(x = sim$final_time-time_grid,
        y = PDRs,
        type = 'l',
        col = 'blue');
}