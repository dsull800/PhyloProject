require("RPANDA")
require("diversitree")
require("castor")
require("phybase")
require("ggtree")
require("treeAGG")
require("ape")

## Not run:
# Generate a random tree with exponentially varying lambda & mu
for(Ntips in c(rep(20,100),rep(100,100),rep(1000,100))){
Ntips = 20
rho = 1 # sampling fraction
time_grid = seq(from=0, to=100, by=0.01)
for(lambdas in list(20+(100/tail(exp(0.1*time_grid),1))*exp(0.1*time_grid),0.2*time_grid+0.5,rep(2,length(time_grid)))){
# lambdas = 2*exp(0.1*time_grid) 44052.93 is the tail of lambdas values
mus = 0*time_grid
sim = castor::generate_random_tree( parameters = list(rarefaction=rho),
                                    max_tips = Ntips/rho,
                                    coalescent = FALSE,
                                    added_rates_times = time_grid,
                                    added_birth_rates_pc = lambdas,added_death_rates_pc = mus)
bigtree = sim[["tree"]]
root_age = castor::get_tree_span(bigtree)[["max_distance"]]
cat(sprintf("Tree has %d tips, spans %g Myr\n",length(bigtree[["tip.label"]]),root_age))
# bigtree[["tip.label"]]=1:length(bigtree[["tip.label"]])
# bigtree[["tip.label"]]=as.integer(bigtree[["tip.label"]])
ape::plot.phylo(bigtree)
title("speciestree")

newickformat=castor::write_tree(bigtree)

# spname=species.name(newickformat)

nodematrix=phybase::read.tree.nodes(newickformat,name=bigtree[["tip.label"]])

#I can set the population sizes by changing the entries of the 5th value of the matrix, 
#use a gamma distributed value or something, maybe.
for(i in 1:nrow(nodematrix[["nodes"]])){
  nodematrix[["nodes"]][i,5]=rgamma(n=1,shape=2,rate=1)
}


nspecies=length(bigtree[["tip.label"]])

# specnames=1:length(bigtree[["tip.label"]])
# specnames=as.integer(bigtree[["tip.label"]])

# specnames=species.name(bigtree)

rootnode=nrow(nodematrix[["nodes"]])

genetreestuff=phybase::sim.coaltree.sp(rootnode,nodematrix[["nodes"]],nspecies,seq=rep(1,nspecies),
                                       name=bigtree[["tip.label"]])

genetreeheight=genetreestuff[["height"]]

realgenetree=genetreestuff[["node"]]

genetreegt=genetreestuff[["gt"]]


thetree=ape::read.tree(text=genetreegt)

ape::plot.phylo(thetree)
title("genetree")

# interimtree=write_tree(thetree)
# thetree=read_tree(interimtree)

# 
# calculate true PDR
lambda_slopes = diff(lambdas)/diff(time_grid);
lambda_slopes = c(lambda_slopes[1],lambda_slopes)
PDRs = lambdas - mus - (lambda_slopes/lambdas)
# Fit PDR on grid
Ngrid = 10
age_grid = seq(0,max(node.depth.edgelength(thetree)),length.out=Ngrid)
# age_grid=seq(0,max(node.depth.edgelength(thetree)),length.out=Ngrid)
fit = castor::fit_hbd_pdr_on_grid(thetree,
                                  age_grid=age_grid,
                                  min_PDR = -Inf,
                                  max_PDR = +Inf,
                                  condition = "crown",
                                  Ntrials = 10,# perform 10 fitting trials
                                  Nthreads = 2,# use two CPUs
                                  max_model_runtime = 1) # limit model evaluation to 1 second
if(!fit[["success"]]){
  cat(sprintf("ERROR: Fitting failed: %s\n",fit[["error"]]))
}else{
  cat(sprintf("Fitting succeeded:\nLoglikelihood=%g\n",fit[["loglikelihood"]]))
  # plot fitted & true PDR
  plot( x = fit[["age_grid"]],
        y = fit[["fitted_PDR"]],
        main = 'Fitted & true PDR',
        xlab = 'age',
        ylab = 'PDR',
        type = 'b',
        col = 'red',
         xlim = c(genetreeheight,0),
        ylim=c(-100,100))
        # xlim = c(-100,100))
  lines(x = seq(from=genetreeheight,to=0,length.out=length(time_grid)),
        y = PDRs,
        type = 'l',
        col = 'blue');
}

# Fit PSR on grid
oldest_age=genetreeheight/2 # only consider recent times when fitting
Ngrid = Ntips
age_grid = seq(from=0,to=oldest_age,length.out=Ngrid)
fit = fit_hbd_psr_on_grid(thetree,
                          oldest_age = oldest_age,
                          age_grid = age_grid,
                          min_PSR = -Inf,
                          max_PSR = +Inf,
                          guess_PSR=1,
                          # fixed_PSR=rep(1,Ngrid),
                          condition = "crown",
                          Ntrials = 10,# perform 10 fitting trials
                          Nthreads = 5,# use two CPUs
                          max_model_runtime = 1) # limit model evaluation to 1 second
if(!fit[["success"]]){
  cat(sprintf("ERROR: Fitting failed: %s\n",fit[["error"]]))
}else{
  cat(sprintf("Fitting succeeded:\nLoglikelihood=%g\n",fit[["loglikelihood"]]))
  # plot fitted PSR
  plot( x = fit[["age_grid"]],
        y = fit[["fitted_PSR"]],
        main = 'Fitted PSR',
        xlab = 'age',
        ylab = 'PSR',
        type = 'b',
        xlim = c(genetreeheight/2,0))
  # plot deterministic LTT of fitted model
  plot( x = fit[["age_grid"]],
        y = fit[["fitted_LTT"]],
        main = 'Fitted dLTT',
        xlab = 'age',
        ylab = 'lineages',
        type = 'b',
        log = 'y',
        xlim = c(genetreeheight/2,0))
}
}#close functions loop
rm(list = ls())
}#close Ntips loop