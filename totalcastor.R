require("castor")
# require("ggplot2")
# require("stats")
require("prospectr")
require("plot.matrix")
require("matrixStats")
require("naniar")


binningstddev <- function(X, bins, bin.size) {
  
  if (is.data.frame(X)) 
    X <- as.matrix(X)
  if (!missing(bins) & !missing(bin.size)) 
    stop("EITHER 'bins' OR 'bin.size' must be specified")
  if (missing(bins) & missing(bin.size)) 
    return(X)
  
  if (is.matrix(X)) 
    p1 <- ncol(X) else p1 <- length(X)
    
    if (missing(bins) & !missing(bin.size)) {
      b <- findInterval(1:p1, seq(1, p1, bin.size))
    } else {
      b <- findInterval(1:p1, seq(1, p1, length.out = bins + 1), rightmost.closed = T)
    }
    
    p2 <- max(b)
    
    if (is.matrix(X)) {
      output <- matrix(0, nrow(X), p2)
      for (i in seq_len(p2)) {
        output[, i] <- rowSds(X[, b == i, drop = F])
      }
      colnames(output) <- colnames(X)[ceiling(tapply(b, b, function(x) mean(which(b == x[1]), na.rm = T)))]  # find colnames
      rownames(output) <- rownames(X)
    } else {
      output <- tapply(X, b, sd)
      names(output) <- names(X)[ceiling(tapply(b, b, function(x) mean(which(b == x[1]), na.rm = T)))]
    }
    
    return(output)
} 





lambdanumber=-1
munumber=-1
Ntipnumber=-1
rho = 1 # sampling fraction
ncols=30
count100=1
count10=1
count1=1
colnamesstuff=c()
for(i in seq(1,ncols)){
  inter_val=toString(i/ncols)
  colnamesstuff=c(colnamesstuff,inter_val)
}
heatmapdata=matrix(0,nrow=3,ncol=ncols)
rownames(heatmapdata)=c("max_val 100","max_val 10","max_val 1")
colnames(heatmapdata)=colnamesstuff
heatmapdatasd=matrix(0,nrow=3,ncol=ncols)
rownames(heatmapdatasd)=c("max_val 100","max_val 10","max_val 1")
colnames(heatmapdatasd)=colnamesstuff

matrix100=matrix(nrow=21*10/3,ncol=ncols)
matrix10=matrix(nrow=21*10/3,ncol=ncols)
matrix1=matrix(nrow=21*10/3,ncol=ncols)

heatmapdatanew=matrix(0,nrow=3,ncol=ncols)
rownames(heatmapdatanew)=c("max_val 100","max_val 10","max_val 1")
colnames(heatmapdatanew)=colnamesstuff
heatmapdatasdnew=matrix(0,nrow=3,ncol=ncols)
rownames(heatmapdatasdnew)=c("max_val 100","max_val 10","max_val 1")
colnames(heatmapdatasdnew)=colnamesstuff

matrix100new=matrix(nrow=21*10/3,ncol=ncols)
matrix10new=matrix(nrow=21*10/3,ncol=ncols)
matrix1new=matrix(nrow=21*10/3,ncol=ncols)

## Rvals are 10,100,1000

for(Ntips in c(rep(20000,21))){
  Ntipnumber=Ntipnumber+1
  if(Ntipnumber%%3==0){
    max_val=100
  }else if(Ntipnumber%%3==1){
    max_val=10
  }else {
    max_val=1
  }
  
  
  time_grid = seq(from=0, to = 1.2*round(log(Ntips/rho)/max_val,digits=1), by=0.01)
  # ln(max_tips)/(lambda-mu)
  lambda1 = exp(0.5*time_grid)
  # for(lambdas in list(((max_val/2)/tail(lambda1,1))*lambda1+max_val/2,(max_val/2*tail(time_grid,1))*time_grid+max_val/2,rep(max_val,length(time_grid)),-(max_val/2/tail(time_grid,1))*time_grid+max_val)){
  #look at constant lambda, then it doesn;t matter which end of the array is which
  for(lambdas in list(((max_val/2)/tail(lambda1,1))*lambda1+max_val/2)){
    A=1.1*lambdas[floor(length(time_grid)/2)]
    sigma=10^-3
    lambdanumber=lambdanumber+1
    # for(mus in list(0*time_grid,A*exp(-(time_grid-time_grid[floor(length(time_grid)/2)])^2/(2*sigma^2)),rep(max_val/3,length(time_grid)))){
    for(mus in list(0*time_grid)){
      munumber=munumber+1
      tryCatch({
        sim = castor::generate_random_tree( parameters = list(rarefaction=rho),
                                            # max_tips = Ntips/rho,
                                            # as_generations = TRUE,
                                            max_time=time_grid[length(time_grid)],
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
        # Ntipnumber was ntips in old simulation
        file=paste(munumber,"_",Ntipnumber,"_",lambdanumber%%4,"_",munumber%%3,".txt",sep="")
        setwd("spectrees")
        file.create(file)
        write_tree(spectree,file)
        setwd("..")
        # 
        #redefine lambdas & mus w.r.t. age_grid
        age_grid = rev(sim$final_time-time_grid)
        lambdas_on_age_grid = rev(lambdas)
        mus_on_age_grid = rev(mus)
        # if extinction is 0 why isn;t sim$final_time=root_age?
        spectreepdrpsr = simulate_deterministic_hbd(LTT0 = length(spectree[["tip.label"]]),
                                                    oldest_age = sim$final_time,
                                                    age0=0,
                                                    age_grid=age_grid,
                                                    rho0 = rho,
                                                    lambda=lambdas_on_age_grid,mu=mus_on_age_grid)
        
        
        #want to simulate multiple genetrees for a given species tree
        genetreenum=-1
        
        for(i in seq(1,10)){
          genetreenum=genetreenum+1
          genetreestuff = generate_gene_tree_msc(spectree,allele_counts = 1,
                                                 population_sizes = 10^8,
                                                 generation_times = 10^-7,
                                                 #runif(nspecies+Nnodes,min = 0.001,max=1) runif(nspecies+Nnodes,min = 10^8,max=10^9)
                                                 ploidy = 1);
          
          gentree=genetreestuff$tree
          
          ape::plot.phylo(gentree)
          title("genetree")
          
          
          # Fit PSR on grid
          #maybe iterate over fineness of grid points?
          Ngrid = 5
          #only fit on species tree height?
          # height=max(get_all_distances_to_root(gentree))
          # height=root_age
          height=sim$final_time
          oldest_age=height # only consider recent times when fitting
          psr_age_grid = seq(from=0,to=oldest_age,length.out=Ngrid)
          fit = fit_hbd_psr_on_grid(gentree,
                                    oldest_age = oldest_age,
                                    age_grid = psr_age_grid,
                                    age0=0,
                                    min_PSR = -50,
                                    max_PSR = +150,
                                    guess_PSR= max_val,
                                    condition = "stem",
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
            # plot deterministic LTT of fitted model, something is weird with the values/plot?
            plot( x = fit[["age_grid"]],
                  y = fit[["fitted_LTT"]],
                  main = 'Fitted dLTT',
                  xlab = 'age',
                  ylab = 'lineages',
                  type = 'b',
                  log = 'y',
                  xlim = c(oldest_age,0))
          }
          
          NGtips = length(gentree$tip.label)
          gene_root_age = castor::get_tree_span(gentree)$max_distance
          gene_LTT = castor::count_lineages_through_time(gentree, max_time=gene_root_age*0.999, Ntimes=floor(sqrt(NGtips)/10), include_slopes=TRUE, regular_grid=FALSE)
          gene_PSR = gene_LTT$relative_slopes
          
          plot(gene_LTT$times, gene_LTT$lineages, type="l", xlab="time", ylab="# clades")
          title("species tree/gene tree LTT")
          
          
          lttcountspec=castor::count_lineages_through_time(spectree,Ntimes=100)
          lines(lttcountspec$times, lttcountspec$lineages, type="l", xlab="time", ylab="# clades")
          
          
          
          ##plot epsilon over time
          lambda_hat_p_prime=approx(x=fit[["age_grid"]],y=fit[["fitted_PSR"]],xout=spectreepdrpsr$ages,method="linear")$y
          lambda_hat_p_prime_new=approx(x=gene_LTT$times,y=gene_PSR,xout=spectreepdrpsr$ages,method="linear")$y
          
          #spectreepdrpsr$PSR is very close to 0, need to use double precision?
          epsilon=(lambda_hat_p_prime-spectreepdrpsr$PSR)/spectreepdrpsr$PSR
          
          NAend=1
          #|NAend==length(epsilon)
          while(is.na(epsilon[NAend])){
            
            NAend=NAend+1
            
          }

          almostrealepsilonvals=epsilon[NAend:length(epsilon)]

          NAend2=1
          # |NAend2==length(almostrealepsilonvals)
          while(!is.na(almostrealepsilonvals[NAend2])){

            NAend2=NAend2+1

          }

          realepsilonvals=almostrealepsilonvals[1:NAend2-1]
          
          if(any_na(realepsilonvals)){
            print(file)
            print("this one has NA realepislonvals")
          }
          

          realepsilonages=spectreepdrpsr$ages[NAend:NAend2-1]/spectreepdrpsr$ages[NAend2-1]
          
          binnedepsilons=binning(realepsilonvals,ncols)
          
          binnedepsilonssd=binningstddev(realepsilonvals,bins=ncols)
          
          
          
          
          epsilonnew=(lambda_hat_p_prime_new-spectreepdrpsr$PSR)/spectreepdrpsr$PSR
          
          NAend=1
          #|NAend==length(epsilon)
          while(is.na(epsilonnew[NAend])){
            
            NAend=NAend+1
            
          }
          
          almostrealepsilonvalsnew=epsilonnew[NAend:length(epsilon)]
          
          NAend2=1
          # |NAend2==length(almostrealepsilonvals)
          while(!is.na(almostrealepsilonvalsnew[NAend2])){
            
            NAend2=NAend2+1
            
          }
          
          realepsilonvalsnew=almostrealepsilonvalsnew[1:NAend2-1]
          
          if(any_na(realepsilonvalsnew)){
            print(file)
            print("this one has NA realepislonvalsnew")
          }
          
          
          realepsilonagesnew=spectreepdrpsr$ages[NAend:NAend2-1]/spectreepdrpsr$ages[NAend2-1]
          
          binnedepsilonsnew=binning(realepsilonvalsnew,ncols)
          
          binnedepsilonssd=binningstddev(realepsilonvalsnew,bins=ncols)
          
          #need to figure out what to divide by to normalize
          if(Ntipnumber%%3==0){
            for(i in 1:length(binnedepsilons)){
              heatmapdata[1,i]=heatmapdata[1,i]+binnedepsilons[i]
              heatmapdatasd[1,i]=heatmapdatasd[1,i]+binnedepsilonssd[i]
              matrix100[count100,i]=binnedepsilons[i]
              
              heatmapdatanew[1,i]=heatmapdatanew[1,i]+binnedepsilonsnew[i]
              heatmapdatasdnew[1,i]=heatmapdatasdnew[1,i]+binnedepsilonssdnew[i]
              matrix100new[count100,i]=binnedepsilonsnew[i]
            }
            count100=count100+1
          }else if(Ntipnumber%%3==1){
            for(i in 1:length(binnedepsilons)){
              heatmapdata[2,i]=heatmapdata[2,i]+binnedepsilons[i]
              heatmapdatasd[2,i]=heatmapdatasd[2,i]+binnedepsilonssd[i]
              matrix10[count10,i]=binnedepsilons[i]
              
              heatmapdatanew[2,i]=heatmapdatanew[2,i]+binnedepsilonsnew[i]
              heatmapdatasdnew[2,i]=heatmapdatasdnew[2,i]+binnedepsilonssdnew[i]
              matrix10new[count10,i]=binnedepsilons[i]
            }
            count10=count10+1
          }else {
            for(i in 1:length(binnedepsilons)){
              heatmapdata[3,i]=heatmapdata[3,i]+binnedepsilons[i]
              heatmapdatasd[3,i]=heatmapdatasd[3,i]+binnedepsilonssd[i]
              matrix1[count1,i]=binnedepsilons[i]
              
              heatmapdatanew[3,i]=heatmapdatanew[3,i]+binnedepsilonsnew[i]
              heatmapdatasdnew[3,i]=heatmapdatasdnew[3,i]+binnedepsilonssdnew[i]
              matrix1new[count1,i]=binnedepsilonsnew[i]
            }
            count1=count1+1
          }
          
          #linear epsilon graph artifact of approx linear interp? when set to constant, epsilon graph is constant, so maybe
          plot(y=epsilon,x=spectreepdrpsr$ages)
          title("epsilon vs. time")
          ### I need to write information to file for each run, I need to index file name by index. Need to include newick strings for gene and species trees, and maybe fitted values for the pdr/psr.
          file=paste(munumber,"_",Ntipnumber,"_",lambdanumber%%4,"_",munumber%%3,"_",genetreenum,".txt",sep="")
          setwd("gentrees")
          file.create(file)
          write_tree(gentree,file)
          
          setwd("../fitpsrs")
          
          file=paste(munumber,"_",Ntipnumber,"_",lambdanumber%%4,"_",munumber%%3,"_",genetreenum,".rds",sep="")
          
          saveRDS(fit,file)
          
          # setwd("../fitpdrs")
          #
          # saveRDS(fitpdr,file)
          
          setwd("../spectreeinfo")
          
          saveRDS(spectreepdrpsr,file)
          
          setwd("../epsilonvals")
          
          saveRDS(binnedepsilons,file)
          
          setwd("../epsilonsd")
          
          saveRDS(binnedepsilonssd,file)
          
          setwd("..")
          
        }#close genetreeloop
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n",file)})
    }#close mus loop
  }#close lambdas loop
}#close Ntips loop

epsilonsd100=colSds(matrix100)
epsilonsd10=colSds(matrix10)
epsilonsd1=colSds(matrix1)

epsilonsdreal=rbind(epsilonsd100,epsilonsd10,epsilonsd1)

plot(epsilonsdreal)
title("sd epsilons")

plot(heatmapdata)
title("average epsilons")

plot(heatmapdatasd)
title("average stddev epsilons")

epsilonsd100new=colSds(matrix100new)
epsilonsd10new=colSds(matrix10new)
epsilonsd1new=colSds(matrix1new)

epsilonsdrealnew=rbind(epsilonsd100new,epsilonsd10new,epsilonsd1new)

plot(epsilonsdrealnew)
title("sd epsilons new")

plot(heatmapdatanew)
title("average epsilons new")

plot(heatmapdatasdnew)
title("average stddev epsilons new")

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="../rundata/plots")














# # calculate true PDR,
# lambda_slopes = diff(lambdas)/diff(time_grid);
# lambda_slopes = c(lambda_slopes[1],lambda_slopes)
# PDRs = lambdas - mus - (lambda_slopes/lambdas)
# # Fit PDR on grid
# Ngrid = 10
# height=max(get_all_distances_to_root(gentree))
# age_grid = seq(0,height,length.out=Ngrid)
# ERROR: Fitting failed: Provided age-grid range (0 - 1.20583) does not cover entire required age range (0 - 1.20583) look at lines 50-51 in fit_hbd_pdr_on_grid
# fitpdr = castor::fit_hbd_pdr_on_grid(gentree,
#                                      age_grid=age_grid,
#                                      min_PDR = -50,
#                                      max_PDR = +150,
#                                      # guess_PDR = tail(PDRs,1),
#                                      condition = "stem",
#                                      Ntrials = 10,# perform 10 fitting trials
#                                      Nthreads = 4,# use two CPUs
#                                      max_model_runtime = 1) # limit model evaluation to 1 second max(1,Ntip/100,000)
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
#         ylim=c(-50,150))
# 
#   lines(x = seq(from=height,to=0,length.out=length(time_grid)),
#         y = PDRs,
#         type = 'l',
#         col = 'blue');
# }