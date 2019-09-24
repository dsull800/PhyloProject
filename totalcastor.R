require("castor")

## Not run:
# Generate a random tree with exponentially varying lambda & mu
for(Ntips in c(rep(20,100),rep(100,100),rep(1000,100))){
# Ntips=20
rho = 1 # sampling fraction
time_grid = seq(from=0, to=100, by=0.01)
for(lambdas in list(((10^-5)/tail(exp(10^-5*time_grid),1))*exp(10^-5*time_grid),10^-6*time_grid,rep(10^-5,length(time_grid)))){
# lambdas=0.2*time_grid+0.5
mus = 0*time_grid
sim = castor::generate_random_tree( parameters = list(rarefaction=rho),
                                    max_tips = Ntips/rho,
                                    # as_generations = TRUE,
                                    coalescent = FALSE,
                                    added_rates_times = time_grid,
                                    added_birth_rates_pc = lambdas,added_death_rates_pc = mus)
bigtree = sim[["tree"]]
root_age = castor::get_tree_span(bigtree)[["max_distance"]]
cat(sprintf("Tree has %d tips, spans %g Myr\n",length(bigtree[["tip.label"]]),root_age))
ape::plot.phylo(bigtree)
title("speciestree")

Nnodes=bigtree$Nnode

nspecies=length(bigtree$tip.label)

# 
genetreestuff = generate_gene_tree_msc(bigtree,allele_counts = 1,
                                       population_sizes = 5*10^8,
                                       generation_times = 0.01,
                                       #runif(nspecies+Nnodes,min = 0.001,max=1) runif(nspecies+Nnodes,min = 10^8,max=10^9)
                                       ploidy = 1);

thetree=genetreestuff$tree

ape::plot.phylo(thetree)
title("genetree")


# calculate true PDR,
lambda_slopes = diff(lambdas)/diff(time_grid);
lambda_slopes = c(lambda_slopes[1],lambda_slopes)
PDRs = lambdas - mus - (lambda_slopes/lambdas)
# Fit PDR on grid
Ngrid = 10
height=max(get_all_distances_to_root(thetree))
age_grid = seq(0,height,length.out=Ngrid)
# ERROR: Fitting failed: Provided age-grid range (0 - 1.20583) does not cover entire required age range (0 - 1.20583) look at lines 50-51 in fit_hbd_pdr_on_grid
fitpdr = castor::fit_hbd_pdr_on_grid(thetree,
                                     age_grid=age_grid,
                                     min_PDR = -20,
                                     max_PDR = +50,
                                     fixed_rholambda0 = rho*tail(lambdas,1),
                                     condition = "crown",
                                     Ntrials = 10,# perform 10 fitting trials
                                     Nthreads = 4,# use two CPUs
                                     max_model_runtime = 1) # limit model evaluation to 1 second
if(!fitpdr[["success"]]){
  cat(sprintf("ERROR: Fitting failed: %s\n",fitpdr[["error"]]))
  stop()
}else{
  cat(sprintf("Fitting succeeded:\nLoglikelihood=%g\n",fitpdr[["loglikelihood"]]))
  # plot fitted & true PDR
  plot( x = fitpdr[["age_grid"]],
        y = fitpdr[["fitted_PDR"]],
        main = 'Fitted & true PDR',
        xlab = 'age',
        ylab = 'PDR',
        type = 'b',
        col = 'red',
        xlim = c(height,0),
        ylim=c(-5,5))
  # xlim = c(-100,100))
  lines(x = seq(from=height,to=0,length.out=length(time_grid)),
        y = PDRs,
        type = 'l',
        col = 'blue');
}

# Fit PSR on grid
oldest_age=height/2 # only consider recent times when fitting
Ngrid = Ntips
age_grid = seq(from=0,to=oldest_age,length.out=Ngrid)
fit = fit_hbd_psr_on_grid(thetree,
                          oldest_age = oldest_age,
                          age_grid = age_grid,
                          min_PSR = -50,
                          max_PSR = +50,
                          guess_PSR=10^-5,
                          condition = "crown",
                          Ntrials = 10,# perform 10 fitting trials
                          Nthreads = 4,# use two CPUs
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
        xlim = c(oldest_age,0))
  # plot deterministic LTT of fitted model
  plot( x = fit[["age_grid"]],
        y = fit[["fitted_LTT"]],
        main = 'Fitted dLTT',
        xlab = 'age',
        ylab = 'lineages',
        type = 'b',
        log = 'y',
        xlim = c(oldest_age,0))
}
lttcount=castor::count_lineages_through_time(thetree,Ntimes=100)
plot(lttcount$times, lttcount$lineages, type="l", xlab="time", ylab="# clades")
  }#close functions loop
  rm(list = ls())
}#close Ntips loop