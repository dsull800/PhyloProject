require("castor")

runnumber=0
## Not run:
# Generate a random tree with exponentially varying lambda & mu
for(Ntips in c(rep(20000,60),rep(50000,60),rep(100000,60))){
rho = 1 # sampling fraction
if(runnumber%%3==0){
  max_val=100
}else if(runnumber%%3==1){
  max_val=10
}else {
  max_val=1
}

time_grid = seq(from=0, to = 2*ceiling(log(Ntips/rho)/max_val), by=0.1)
# ln(max_tips)/(lambda-mu)
lambda1 = exp(.1*time_grid)

for(lambdas in list(((max_val)/tail(lambda1,1))*lambda1,(max_val/tail(time_grid,1))*time_grid,rep(max_val,length(time_grid)))){
tryCatch({
# print(tail(lambda1,1))
# lambdas=((100)/tail(lambda1,1))*lambda1
mus = 0*time_grid
sim = castor::generate_random_tree( parameters = list(rarefaction=rho),
                                    max_tips = Ntips/rho,
                                    # as_generations = TRUE,
                                    coalescent = TRUE,
                                    added_rates_times = time_grid,
                                    added_birth_rates_pc = lambdas,added_death_rates_pc = mus)



###ALL UNITS ARE IN MEGAYEARS
###ALL UNITS ARE IN MEGAYEARS

spectree = sim[["tree"]]
root_age = castor::get_tree_span(spectree)[["max_distance"]]
cat(sprintf("Tree has %d tips, spans %g Myr\n",length(spectree[["tip.label"]]),root_age))
ape::plot.phylo(spectree)
title("speciestree")

Nnodes=spectree$Nnode

nspecies=length(spectree$tip.label)

# 
genetreestuff = generate_gene_tree_msc(spectree,allele_counts = 1,
                                       population_sizes = 10^8,
                                       generation_times = 10^5,
                                       #runif(nspecies+Nnodes,min = 0.001,max=1) runif(nspecies+Nnodes,min = 10^8,max=10^9)
                                       ploidy = 1);

gentree=genetreestuff$tree

ape::plot.phylo(gentree)
title("genetree")


# calculate true PDR,
lambda_slopes = diff(lambdas)/diff(time_grid);
lambda_slopes = c(lambda_slopes[1],lambda_slopes)
PDRs = lambdas - mus - (lambda_slopes/lambdas)
# Fit PDR on grid
Ngrid = 10
height=max(get_all_distances_to_root(gentree))
age_grid = seq(0,height,length.out=Ngrid)
# ERROR: Fitting failed: Provided age-grid range (0 - 1.20583) does not cover entire required age range (0 - 1.20583) look at lines 50-51 in fit_hbd_pdr_on_grid
fitpdr = castor::fit_hbd_pdr_on_grid(gentree,
                                     age_grid=age_grid,
                                     min_PDR = -50,
                                     max_PDR = +200,
                                     fixed_rholambda0 = rho*tail(lambdas,1),
                                     # guess_PDR = tail(PDRs,1),
                                     condition = "crown",
                                     Ntrials = 10,# perform 10 fitting trials
                                     Nthreads = 12,# use two CPUs
                                     max_model_runtime = 10) # limit model evaluation to 1 second
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
        ylim=c(-5,100))
 
  lines(x = seq(from=height,to=0,length.out=length(time_grid)),
        y = PDRs,
        type = 'l',
        col = 'blue');
}

# Fit PSR on grid
oldest_age=height/2 # only consider recent times when fitting
Ngrid = Ntips
age_grid = seq(from=0,to=oldest_age,length.out=Ngrid)
fit = fit_hbd_psr_on_grid(gentree,
                          oldest_age = oldest_age,
                          age_grid = age_grid,
                          min_PSR = -50,
                          max_PSR = +200,
                          guess_PSR= tail(lambdas,1),
                          condition = "crown",
                          Ntrials = 10,# perform 10 fitting trials
                          Nthreads = 12,# use two CPUs
                          max_model_runtime = 10) # limit model evaluation to 1 second
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
lttcountgen=castor::count_lineages_through_time(gentree,Ntimes=100)
plot(lttcountgen$times, lttcountgen$lineages, type="l", xlab="time", ylab="# clades")

lttcountspec=castor::count_lineages_through_time(spectree,Ntimes=100)
plot(lttcountspec$times, lttcountspec$lineages, type="l", xlab="time", ylab="# clades")
title()

funcnumber=runnumber%%3
runnumber=runnumber+1
### I need to write information to file for each run, I need to index file name by index. Need to include newick strings for gene and species trees, and maybe fitted values for the pdr/psr. 
setwd("gentrees")
file=paste(runnumber,"_",Ntips,"_",funcnumber,".txt",sep="")
file.create(file)
write_tree(gentree,file)

setwd("../spectrees")
file.create(file)
write_tree(spectree,file)

setwd("..")

file=paste(runnumber,"_",Ntips,"_",funcnumber,".rds",sep="")

setwd("fitpsrs")

saveRDS(fit,file)

setwd("../fitpdrs")

saveRDS(fitpdr,file)

setwd("..")
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }#close functions loop
}#close Ntips loop







































# # Fit PDR on grid for the species tree
# Ngrid = 10
# height=max(get_all_distances_to_root(spectree))
# age_grid = seq(0,height,length.out=Ngrid)
# # ERROR: Fitting failed: Provided age-grid range (0 - 1.20583) does not cover entire required age range (0 - 1.20583) look at lines 50-51 in fit_hbd_pdr_on_grid
# fitpdr = castor::fit_hbd_pdr_on_grid(spectree,
#                                      age_grid=age_grid,
#                                      min_PDR = -50,
#                                      max_PDR = +200,
#                                      fixed_rholambda0 = rho*tail(lambdas,1),
#                                      # guess_PDR = tail(PDRs,1),
#                                      condition = "crown",
#                                      Ntrials = 10,# perform 10 fitting trials
#                                      Nthreads = 12,# use two CPUs
#                                      max_model_runtime = 1) # limit model evaluation to 1 second
# if(!fitpdr[["success"]]){
#   cat(sprintf("ERROR: Fitting failed: %s\n",fitpdr[["error"]]))
#   stop()
# }else{
#   cat(sprintf("Fitting succeeded:\nLoglikelihood=%g\n",fitpdr[["loglikelihood"]]))
#   # plot fitted & true PDR
#   plot( x = fitpdr[["age_grid"]],
#         y = fitpdr[["fitted_PDR"]],
#         main = 'Fitted & true PDR',
#         xlab = 'age',
#         ylab = 'PDR',
#         type = 'b',
#         col = 'red',
#         xlim = c(height,0),
#         ylim=c(-5,100))
#   
#   lines(x = seq(from=height,to=0,length.out=length(time_grid)),
#         y = PDRs,
#         type = 'l',
#         col = 'blue');
# }