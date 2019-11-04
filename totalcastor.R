require("castor")
# require("ggplot2")
# require("stats")
require("prospectr")
require("plot.matrix")
require("matrixStats")
require("naniar")


lambdanumber=-1
munumber=-1
Ntipnumber=-1
rho = 1 # sampling fraction
count100=1
count10=1
count1=1
oldest_age_sim=1
age_grid_fineness=.01
ncols=oldest_age_sim/age_grid_fineness
colnamesstuff=c()
for(i in seq(1,ncols)){
  inter_val=toString(i/ncols)
  colnamesstuff=c(colnamesstuff,inter_val)
}
heatmapdata=matrix(0,nrow=3,ncol=ncols)
rownames(heatmapdata)=c("max_val 100","max_val 10","max_val 1")
colnames(heatmapdata)=colnamesstuff

matrix100=matrix(nrow=21*10/3,ncol=ncols)
matrix10=matrix(nrow=21*10/3,ncol=ncols)
matrix1=matrix(nrow=21*10/3,ncol=ncols)

heatmapdatanew=matrix(0,nrow=3,ncol=ncols)
rownames(heatmapdatanew)=c("max_val 100","max_val 10","max_val 1")
colnames(heatmapdatanew)=colnamesstuff

matrix100new=matrix(nrow=21*10/3,ncol=ncols)
matrix10new=matrix(nrow=21*10/3,ncol=ncols)
matrix1new=matrix(nrow=21*10/3,ncol=ncols)

## Rvals are 10,100,1000
#21-11 for ntipnumber
for(Ntips in c(rep(20000,21))){
  Ntipnumber=Ntipnumber+1
  if(Ntipnumber%%3==0){
    max_val=100
  }else if(Ntipnumber%%3==1){
    max_val=10
  }else {
    max_val=1
  }
  
  age_grid_sim = seq(from=0, to = oldest_age_sim, by=age_grid_fineness)
  # ln(max_tips)/(lambda-mu)
  lambda1 = exp(0.5*age_grid_sim)
  # for(lambdas in list(((max_val/2)/tail(lambda1,1))*lambda1+max_val/2,(max_val/2*tail(age_grid_sim,1))*age_grid_sim+max_val/2,rep(max_val,length(age_grid_sim)),-(max_val/2/tail(age_grid_sim,1))*age_grid_sim+max_val)){
  #look at constant lambda, then it doesn;t matter which end of the array is which
  for(lambdas in list(((max_val/2)/tail(lambda1,1))*lambda1+max_val/2)){
    A=1.1*lambdas[floor(length(age_grid_sim)/2)]
    sigma=10^-3
    lambdanumber=lambdanumber+1
    # for(mus in list(0*age_grid_sim,A*exp(-(age_grid_sim-age_grid_sim[floor(length(age_grid_sim)/2)])^2/(2*sigma^2)),rep(max_val/3,length(age_grid_sim)))){
    for(mus in list(A*exp(-(age_grid_sim-age_grid_sim[floor(length(age_grid_sim)/2)])^2/(2*sigma^2)))){
      munumber=munumber+1
      tryCatch({
        sim= castor::generate_tree_hbd_reverse(Ntips=Ntips, age_grid=age_grid_sim, lambda=lambdas,mu=mus,crown_age=oldest_age_sim,rho=rho)
        

        ###ALL UNITS ARE IN MEGAYEARS
        ###ALL UNITS ARE IN MEGAYEARS
        
        spectree = sim$trees[[1]]
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
      
        #redefine lambdas & mus w.r.t. age_grid
        lambdas_on_age_grid = lambdas
        mus_on_age_grid = mus
        # if extinction is 0 why isn;t sim$final_time=root_age? Could I just calculate this once and then store it?
        spectreepdrpsr = simulate_deterministic_hbd(LTT0 = length(spectree[["tip.label"]]),
                                                    # oldest_age = oldest_age_sim,
                                                    # age0=0,
                                                    age_grid=age_grid_sim,
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
          
          Ngrid = 5
          
          psr_age_grid = seq(from=0,to=oldest_age,length.out=Ngrid)
          fit = fit_hbd_psr_on_grid(gentree,
                                    oldest_age = oldest_age_sim,
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
          # gene_root_age = castor::get_tree_span(gentree)$max_distance
          #for max_time use oldest age_sim because all we need to fit is time that species trees exists
          gene_LTT = castor::count_lineages_through_time(gentree, max_time=oldest_age_sim, Ntimes=len(oldest_age_sim), include_slopes=TRUE, regular_grid=TRUE)
          gene_PSR = gene_LTT$relative_slopes
          
          plot(gene_LTT$times, gene_LTT$lineages, type="l", xlab="time", ylab="# clades")
          title("species tree/gene tree LTT")
          
          
          lttcountspec=castor::count_lineages_through_time(spectree,Ntimes=100)
          lines(lttcountspec$times, lttcountspec$lineages, type="l", xlab="time", ylab="# clades")
          
          # lttcountspec$relative_slopes[n] = PSR at time lttcountspec$times[n] and thus at age root_age-lttcountspec$times[n]
          # --> lttcountspec$relative_slopes[] is synchronized with ages[] = root_age - lttcountspec$times[]
          
          
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

          binnedepsilons=realepsilonvals
          

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
          
          binnedepsilonsnew=realepsilonvalsnew
          
          #need to figure out what to divide by to normalize
          if(Ntipnumber%%3==0){
            for(i in 1:length(binnedepsilons)){
              heatmapdata[1,i]=heatmapdata[1,i]+binnedepsilons[i]
              
              matrix100[count100,i]=binnedepsilons[i]
            
            }
            
            for(i in 1:length(binnedepsilonsnew)){
              heatmapdatanew[1,i]=heatmapdatanew[1,i]+binnedepsilonsnew[i]
              
              matrix100new[count100,i]=binnedepsilonsnew[i]
            }
            count100=count100+1
          }else if(Ntipnumber%%3==1){
            for(i in 1:length(binnedepsilons)){
              heatmapdata[2,i]=heatmapdata[2,i]+binnedepsilons[i]
              
              matrix10[count10,i]=binnedepsilons[i]
              
            
            }
            
            for(i in 1:length(binnedepsilonsnew)){
              heatmapdatanew[2,i]=heatmapdatanew[2,i]+binnedepsilonsnew[i]
              
              matrix10new[count10,i]=binnedepsilons[i]
            }
            count10=count10+1
          }else {
            for(i in 1:length(binnedepsilons)){
              heatmapdata[3,i]=heatmapdata[3,i]+binnedepsilons[i]
              
              matrix1[count1,i]=binnedepsilons[i]
              

            }
            for(i in 1:length(binnedepsilonsnew)){
              heatmapdatanew[3,i]=heatmapdatanew[3,i]+binnedepsilonsnew[i]
              
              matrix1new[count1,i]=binnedepsilonsnew[i]
            }
            
            count1=count1+1
          }
          file=paste(munumber,"_",Ntipnumber,"_",lambdanumber%%4,"_",munumber%%3,"_",genetreenum,".pdf",sep="")
          
          #linear epsilon graph artifact of approx linear interp? when set to constant, epsilon graph is constant, so maybe
          ##changeworking DIRECTORY
          pdf(file=file, width=5, height=5)
          plot(y=realepsilonvals,x=spectreepdrpsr$ages)
          title("epsilon vs. time")
          
          
          plot(y=realepsilonvalsnew,x=spectreepdrpsr$ages)
          title("realepsilon vs. time")
          
          invisible(dev.off());
          
          ### I need to write information to file for each run, I need to index file name by index. Need to include newick strings for gene and species trees, and maybe fitted values for the pdr/psr.
          file=paste(munumber,"_",Ntipnumber,"_",lambdanumber%%4,"_",munumber%%3,"_",genetreenum,".txt",sep="")
          setwd("gentrees")
          file.create(file)
          write_tree(gentree,file)
          
          setwd("../fitpsrs")
          
          file=paste(munumber,"_",Ntipnumber,"_",lambdanumber%%4,"_",munumber%%3,"_",genetreenum,".rds",sep="")
          
          saveRDS(fit,file)
          
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

par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(epsilonsdreal)
title("sd epsilons")

par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(heatmapdata)
title("average epsilons")


epsilonsd100new=colSds(matrix100new)
epsilonsd10new=colSds(matrix10new)
epsilonsd1new=colSds(matrix1new)

epsilonsdrealnew=rbind(epsilonsd100new,epsilonsd10new,epsilonsd1new)

par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(epsilonsdrealnew)
title("sd epsilons new")

par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(heatmapdatanew)
title("average epsilons new")
