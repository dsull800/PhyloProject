##Something is wrong with the gene_PSR, need to account for crown_age=20 when making plots and doing other calculations
#plots go from present at the left to past at the right

#If very large trees are deterministic, should I only simulate one species tree for given R value?

require("castor")
require("prospectr")
require("plot.matrix")
require("matrixStats")
require("naniar")
overallcount=0
#loop through functions for different scenarios
for(age2lambda in c(function(ages) rep(1,length(ages)),function(ages) rev(ages/max(ages)))){
  for(age2mu in c(function(ages) rep(0,length(ages)), function(ages) 0.1 + 1*exp(-(ages-ages[floor(length(ages)/2)])^2/(2*0.5^2)))){
    overallcount=overallcount+1
    #set working directory depedning on loop variable
    setwd(paste("/Users/danielsullivan/desktop/phylobashstuff/PhyloProject-masterv",toString(overallcount),sep=""))
    #set some loop variables for matrices and variables for grids
    count100=1
    count10=1
    count1=1
    countspec=1
    oldest_age_sim=1000
    #set fineness of age grid
    age_grid_fineness=.1
    #make column names for matrices
    ncols=oldest_age_sim/age_grid_fineness+1
    colnamesstuff=c()
    for(i in seq(1,ncols)){
      inter_val=toString(i/ncols)
      colnamesstuff=c(colnamesstuff,inter_val)
    }
    #create matrices for storing info
    heatmapdata=matrix(0,nrow=3,ncol=ncols)
    rownames(heatmapdata)=c("max_val 100","max_val 10","max_val 1")
    colnames(heatmapdata)=colnamesstuff
    
    spectreematrix=matrix(0,nrow=21,ncol=ncols)
    
    matrix100=matrix(0,nrow=21*10/3,ncol=ncols)
    matrix10=matrix(0,nrow=21*10/3,ncol=ncols)
    matrix1=matrix(0,nrow=21*10/3,ncol=ncols)
    
    heatmapdatanew=matrix(0,nrow=3,ncol=ncols)
    rownames(heatmapdatanew)=c("max_val 100","max_val 10","max_val 1")
    colnames(heatmapdatanew)=colnamesstuff
    
    matrix100new=matrix(0,nrow=21*10/3,ncol=ncols)
    matrix10new=matrix(0,nrow=21*10/3,ncol=ncols)
    matrix1new=matrix(0,nrow=21*10/3,ncol=ncols)
    #initialize loop variable
    Ntipnumber=-1
    
    for(Ntips in c(rep(100000,21))){
      Ntipnumber=Ntipnumber+1
      if(Ntipnumber%%3==0){
     
        max_val=10
      }else if(Ntipnumber%%3==1){
       
        max_val=1
      }else {
   
        max_val=10^-1
      }
      
      #set age grid and rho
      age_grid_sim = seq(from=0, to = oldest_age_sim, by=age_grid_fineness)
      rho			= .5
      #simulate trees
      for(lambdas in list(age2lambda(age_grid_sim))){
        for(mus in list(age2mu(age_grid_sim))){
          # tryCatch({
          
          findcrownage = simulate_deterministic_hbd(LTT0 = Ntips,
                                                    oldest_age = oldest_age_sim,
                                                    age0=0,
                                                    age_grid=age_grid_sim,
                                                    rho0 = rho,
                                                    lambda=lambdas,mu=mus)
          
          crown_age=max(findcrownage$ages[findcrownage$LTT>=1])
          
          lineagecountgrid=seq(from=0,to=crown_age,length.out=1000)
          
          
          sim= castor::generate_tree_hbd_reverse(Ntips=Ntips, age_grid=age_grid_sim, lambda=lambdas,mu=mus,crown_age=crown_age,rho=rho)
          
          
          ###ALL UNITS ARE IN MEGAYEARS
          ###ALL UNITS ARE IN MEGAYEARS
          
          spectree = sim$trees[[1]]
          root_age = castor::get_tree_span(spectree)[["max_distance"]]
          cat(sprintf("Tree has %d tips, spans %g Myr\n",length(spectree[["tip.label"]]),root_age))
          
          
          file=paste(Ntipnumber,"_",".txt",sep="")
          setwd("spectrees")
          file.create(file)
          write_tree(spectree,file)
          setwd("..")
          
          # if extinction is 0 why isn;t sim$final_time=root_age? Could I just calculate this once and then store it?
          spectreepdrpsr = simulate_deterministic_hbd(LTT0 = length(spectree[["tip.label"]]),
                                                      oldest_age = crown_age,
                                                      age0=0,
                                                      age_grid=age_grid_sim,
                                                      rho0 = rho,
                                                      lambda=lambdas,mu=mus)
          
          lttcountspec=castor::count_lineages_through_time(spectree,times=lineagecountgrid,include_slopes = TRUE)
          
          #maybe should make a plot of the relative error between the PSRs of spectree and sim deterministic hbd 
          real_lambda_hat=approx(x=spectreepdrpsr$ages,y=spectreepdrpsr$PSR,xout=lineagecountgrid,method="linear")$y
          
          lambda_hat_spectree=rev(lttcountspec$relative_slopes)
          
          epsilonspectree=(lambda_hat_spectree-real_lambda_hat)/real_lambda_hat
          
          setwd("lambdaplots")
          
          file=paste(Ntipnumber,".pdf",sep="")
          pdf(file=file, width=5, height=5)
          plot(y=real_lambda_hat,x=lineagecountgrid,ylim=c(min(real_lambda_hat),max(real_lambda_hat)))
          plot(y=lambda_hat_spectree,x=lineagecountgrid,ylim=c(min(real_lambda_hat),max(real_lambda_hat)))
          invisible(dev.off());
          
          file=paste(Ntipnumber,".rds",sep="")
          saveRDS(c(real_lambda_hat,lambda_hat_spectree),file)
          
          setwd("..")
          
          for(i in 1:length(epsilonspectree)){
            spectreematrix[countspec,i]=spectreematrix[countspec,i]+epsilonspectree[i]
          }
          countspec=countspec+1
          
          #want to simulate multiple genetrees for a given species tree
          genetreenum=-1
          
          for(i in seq(1,10)){
            genetreenum=genetreenum+1
            genetreestuff = generate_gene_tree_msc(spectree,allele_counts = 1,
                                                   population_sizes = 10^8/max_val,
                                                   generation_times = 10^-7,
                                                   ploidy = 1);
            
            gentree=genetreestuff$tree
            
            # Fit PSR on grid
            
            Ngrid = 10
            
            psr_age_grid = seq(from=0,to=crown_age,length.out=Ngrid)
            fit = fit_hbd_psr_on_grid(gentree,
                                      oldest_age = crown_age,
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
                    type = 'b'
                    # xlim = c(oldest_age_sim,0)
              )
              # plot deterministic LTT of fitted model, something is weird with the values/plot?
              plot( x = fit[["age_grid"]],
                    y = fit[["fitted_LTT"]],
                    main = 'Fitted dLTT',
                    xlab = 'age',
                    ylab = 'lineages',
                    type = 'b',
                    log = 'y'
                    # xlim = c(oldest_age_sim,0)
              )
            }
            
            root_age = castor::get_tree_span(spectree)[["max_distance"]]
            root_age_gene_tree=castor::get_tree_span(gentree)[["max_distance"]]
            distancebetween=root_age_gene_tree-root_age
            #need to increase fineness of times so that finite difference errors are mitigated, actually need to remove dependence on age_grid_sim
            gene_LTT = castor::count_lineages_through_time(gentree, 
                                                           # max_time=oldest_age_sim, Ntimes=length(oldest_age_sim), 
                                                           times=lineagecountgrid+distancebetween,
                                                           include_slopes=TRUE, 
                                                           regular_grid=TRUE)
            gene_PSR = gene_LTT$relative_slopes
            
            #this plot goes from past to present, reverse of what is standard in the rest of the code
            plot(gene_LTT$times-distancebetween, gene_LTT$lineages, type="l", xlab="time", ylab="# clades",col="red",ylim =c(0,101000))
            lines(lttcountspec$times, lttcountspec$lineages, type="l", xlab="time", ylab="# clades",ylim=c(0,101000))
            title("species tree/gene tree LTT")
            
            # lttcountspec$relative_slopes[n] = PSR at time lttcountspec$times[n] and thus at age root_age-lttcountspec$times[n]
            # --> lttcountspec$relative_slopes[] is synchronized with ages[] = root_age - lttcountspec$times[]
            
            
            ##plot epsilon over time
            #xout shouldn;t depend on age_grid_sim, should be less
            lambda_hat_p_prime=approx(x=fit[["age_grid"]],y=fit[["fitted_PSR"]],xout=lineagecountgrid,method="linear")$y
            lambda_hat_p_prime_new=rev(gene_PSR)
            
            #spectreepdrpsr$PSR is very close to 0, need to use double precision?
            epsilon=(lambda_hat_p_prime-real_lambda_hat)/real_lambda_hat
            
            epsilonnew=(lambda_hat_p_prime_new-real_lambda_hat)/real_lambda_hat
            
            
            #need to figure out what to divide by to normalize
            if(Ntipnumber%%3==0){
              for(i in 1:length(epsilon)){
                heatmapdata[1,i]=heatmapdata[1,i]+epsilon[i]
                
                matrix100[count100,i]=epsilon[i]
                
              }
              
              for(i in 1:length(epsilonnew)){
                heatmapdatanew[1,i]=heatmapdatanew[1,i]+epsilonnew[i]
                
                matrix100new[count100,i]=epsilonnew[i]
              }
              count100=count100+1

              
            }else if(Ntipnumber%%3==1){
              for(i in 1:length(epsilon)){
                heatmapdata[2,i]=heatmapdata[2,i]+epsilon[i]
                
                matrix10[count10,i]=epsilon[i]
                
                
              }
              
              for(i in 1:length(epsilonnew)){
                heatmapdatanew[2,i]=heatmapdatanew[2,i]+epsilonnew[i]
                
                matrix10new[count10,i]=epsilon[i]
              }
              count10=count10+1
              
              
              
            }else {
              for(i in 1:length(epsilon)){
                heatmapdata[3,i]=heatmapdata[3,i]+epsilon[i]
                
                matrix1[count1,i]=epsilon[i]
                
                
              }
              for(i in 1:length(epsilonnew)){
                heatmapdatanew[3,i]=heatmapdatanew[3,i]+epsilonnew[i]
                
                matrix1new[count1,i]=epsilonnew[i]
              }
              
              count1=count1+1
              
              
            }
            file=paste(Ntipnumber,"_",genetreenum,".pdf",sep="")
            
            #linear epsilon graph artifact of approx linear interp? when set to constant, epsilon graph is constant, so maybe
            ##changeworking DIRECTORY
            setwd("storedplots")
            pdf(file=file, width=5, height=5)
            plot(y=epsilon,x=lineagecountgrid)
            title("epsilon vs. time")
            
            
            plot(y=epsilonnew,x=lineagecountgrid)
            title("realepsilon vs. time")
            
            plot(y=lambda_hat_p_prime_new,x=lineagecountgrid)
            title("GENE PSR")
            
            invisible(dev.off());
            
            ### I need to write information to file for each run, I need to index file name by index. Need to include newick strings for gene and species trees, and maybe fitted values for the pdr/psr.
            file=paste(Ntipnumber,"_",genetreenum,".txt",sep="")
            setwd("../gentrees")
            file.create(file)
            write_tree(gentree,file)
            
            setwd("../fitpsrs")
            
            file=paste(Ntipnumber,"_",genetreenum,".rds",sep="")
            
            saveRDS(fit,file)
            
            setwd("../spectreeinfo")
            
            saveRDS(spectreepdrpsr,file)
            
            setwd("../epsilonvals")
            
            saveRDS(epsilon,file)
            
            setwd("..")
            
          }#close genetreeloop
          # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n",file)})
        }#close mus loop
      }#close lambdas loop
    }#close Ntips loop
    
    epsilonsd100=colSds(matrix100)
    epsilonsd10=colSds(matrix10)
    epsilonsd1=colSds(matrix1)
    
    epsilonsdreal=rbind(epsilonsd100,epsilonsd10,epsilonsd1)
    
    epsilonsd100new=colSds(matrix100new)
    epsilonsd10new=colSds(matrix10new)
    epsilonsd1new=colSds(matrix1new)
    
    epsilonsdrealnew=rbind(epsilonsd100new,epsilonsd10new,epsilonsd1new)
    
    setwd("matrixplots")
    file=paste(Ntipnumber,"_",genetreenum,".pdf",sep="")
    pdf(file=file, width=5, height=5)
    
    par(mar=c(5.1, 4.1, 4.1, 4.1))
    plot(epsilonsdreal,border=NA,col=hcl.colors(10, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))
    # title("sd epsilons")
    
    par(mar=c(5.1, 4.1, 4.1, 4.1))
    plot(heatmapdata/(21*10/3),border=NA,col=hcl.colors(10, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))
    # title("average epsilons")
    
    par(mar=c(5.1, 4.1, 4.1, 4.1))
    plot(epsilonsdrealnew,border=NA,col=hcl.colors(10, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))
    # title("sd epsilons new")
    
    par(mar=c(5.1, 5.1, 5.1, 5.1))
    plot(heatmapdatanew/(21*10/3),border=NA,col=hcl.colors(10, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))
    #title("average epsilons new")
    
    par(mar=c(5.1, 5.1, 5.1, 5.1))
    plot(spectreematrix,border=NA,col=hcl.colors(50, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))

    invisible(dev.off());
    
    save.image()
    
    setwd("..")
    
  }
}