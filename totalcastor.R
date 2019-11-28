require("castor")
require("prospectr")
require("plot.matrix")
require("matrixStats")
require("naniar")
overallcount=0
#loop through functions for different scenarios
for(age2lambda in c(function(ages) rep(1,length(ages)))){
  for(age2mu in c(function(ages) rep(0,length(ages)), function(ages) 0.1 + 1*exp(-(ages-ages[floor(length(ages)/2)])^2/(2*0.5^2)))){
    overallcount=overallcount+1
    #set working directory depedning on loop variable
    setwd(paste("/Users/danielsullivan/desktop/phylobashstuff/PhyloProject-masterv",toString(overallcount),sep=""))
    #set some loop variables for matrices and variables for grids
    oldest_age_sim=1000
    lineagecount=100
    #set fineness of age grid
    age_grid_fineness=.1
    #make column names for matrices
    rho			= .5
    #more varibales for column names
    ncols=lineagecount
    colnamesstuff=c()
    for(i in seq(1,ncols)){
      inter_val=toString(i/ncols)
      colnamesstuff=c(colnamesstuff,inter_val)
    }
    # vector of R values
    Rvec=c(10^-1,1,10)
    #number of spectrees to generate
    numberofspec=21
    #number of gene trees to generate for each species tree
    numberofgen=10
    
    spectreematrix=matrix(0,nrow=numberofspec,ncol=ncols)
    
    heatmapdatanew=matrix(0,nrow=length(Rvec),ncol=ncols)
    rownames(heatmapdatanew)=Rvec
    colnames(heatmapdatanew)=colnamesstuff
    
    heatmapdatanewsds=matrix(0,nrow=length(Rvec),ncol=ncols)
    rownames(heatmapdatanew)=Rvec
    colnames(heatmapdatanew)=colnamesstuff

    
    for(R in Rvec){
      Rmatrix=matrix(0,nrow=numberofspec*numberofgen,ncol=lineagecount)
      #initialize loop variable
      Ntipnumber=-1
    for(Ntips in c(rep(100000,numberofspec))){
      Ntipnumber=Ntipnumber+1
      
      #set age grid and rho
      age_grid_sim = seq(from=0, to = oldest_age_sim, by=age_grid_fineness)
     
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
          
          lineagecountgrid=seq(from=0,to=crown_age,length.out=lineagecount)
          
          
          sim= castor::generate_tree_hbd_reverse(Ntips=Ntips, age_grid=age_grid_sim, lambda=lambdas,mu=mus,crown_age=crown_age,rho=rho)
          
          
          ###ALL UNITS ARE IN MEGAYEARS
          ###ALL UNITS ARE IN MEGAYEARS
          
          spectree = sim$trees[[1]]
          root_age = castor::get_tree_span(spectree)[["max_distance"]]
          cat(sprintf("Tree has %d tips, spans %g Myr\n",length(spectree[["tip.label"]]),root_age))
          
          
          file=paste(R,"_",Ntipnumber,"_",".txt",sep="")
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
          
          file=paste(R,"_",Ntipnumber,".pdf",sep="")
          pdf(file=file, width=5, height=5)
          plot(y=real_lambda_hat,x=lineagecountgrid,ylim=c(min(real_lambda_hat),max(real_lambda_hat)))
          plot(y=lambda_hat_spectree,x=lineagecountgrid,ylim=c(min(real_lambda_hat),max(real_lambda_hat)))
          invisible(dev.off());
          
          file=paste(R,"_",Ntipnumber,".rds",sep="")
          saveRDS(c(real_lambda_hat,lambda_hat_spectree),file)
          
          setwd("..")
          
          for(i in 1:length(epsilonspectree)){
            spectreematrix[Ntipnumber,i]=spectreematrix[Ntipnumber,i]+epsilonspectree[i]
          }
          
          #want to simulate multiple genetrees for a given species tree
          genetreenum=-1
          
          gentime=10^-7
          popsize=age2lambda(0)/(R*gentime)
          for(i in seq(1,numberofspec)){
            genetreenum=genetreenum+1
            genetreestuff = generate_gene_tree_msc(spectree,allele_counts = 1,
                                                   population_sizes = popsize,
                                                   generation_times = gentime,
                                                   ploidy = 1);
            
            gentree=genetreestuff$tree

            #compute distance between to get values from gene tree for comparisons
            root_age = castor::get_tree_span(spectree)[["max_distance"]]
            root_age_gene_tree=castor::get_tree_span(gentree)[["max_distance"]]
            distancebetween=root_age_gene_tree-root_age
            
            #get genePSR as function of time
            gene_LTT = castor::count_lineages_through_time(gentree, 
                                                           # max_time=oldest_age_sim, Ntimes=length(oldest_age_sim), 
                                                           times=lineagecountgrid+distancebetween,
                                                           include_slopes=TRUE, 
                                                           regular_grid=TRUE)
            #make genePSR as function of age
            gene_PSR = gene_LTT$relative_slopes
            
            #this plot goes from past to present, reverse of what is standard in the rest of the code
            plot(gene_LTT$times-distancebetween, gene_LTT$lineages, type="l", xlab="time", ylab="# clades",col="red",ylim =c(0,101000))
            lines(lttcountspec$times, lttcountspec$lineages, type="l", xlab="time", ylab="# clades",ylim=c(0,101000))
            title("species tree/gene tree LTT")
            
            #get lambdahatpprime values for later
            
            lambda_hat_p_prime_new=rev(gene_PSR)
            
            #compute epsilon values 
            
            epsilonnew=(lambda_hat_p_prime_new-real_lambda_hat)/real_lambda_hat

            # Rmatrix stuff
            
            for(i in seq(1,length(epsilonnew))){
              Rmatrix[(Ntipnumber+1)*10-9+genetreenum,i]=epsilonnew[i]
            }
            
            file=paste(Ntipnumber,"_",genetreenum,".pdf",sep="")
            
            # save plots of epsilon vs age and gene PSR
            setwd("storedplots")
            pdf(file=file, width=5, height=5)
            
            
            plot(y=epsilonnew,x=lineagecountgrid)
            title("realepsilon vs. time")
            
            plot(y=lambda_hat_p_prime_new,x=lineagecountgrid)
            title("GENE PSR")
            
            #this plot goes from past to present, reverse of what is standard in the rest of the code
            plot(gene_LTT$times-distancebetween, gene_LTT$lineages, type="l", xlab="time", ylab="# clades",col="red",ylim =c(0,101000))
            lines(lttcountspec$times, lttcountspec$lineages, type="l", xlab="time", ylab="# clades",ylim=c(0,101000))
            title("species tree/gene tree LTT")
            invisible(dev.off());
            
            # save information to files in certain directories
            file=paste(Ntipnumber,"_",genetreenum,".txt",sep="")
            setwd("../gentrees")
            file.create(file)
            write_tree(gentree,file)

            setwd("../spectreeinfo")
            
            saveRDS(spectreepdrpsr,file)

            setwd("..")
            
          }#close genetreeloop
          # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n",file)})
        }#close mus loop
      }#close lambdas loop
    }#close Ntips loop
      whatisR=Rvec==R
      for(i in seq(1:lineagecount)){
      heatmapdatanew[which(whatisR),i]=colMeans(Rmatrix)[i]
      
      heatmapdatanewsds[which(whatisR),i]=colSds(Rmatrix)[i]
      }
    } #close r loop
    
    #create and plot matrices
    setwd("matrixplots")
    file=paste(Ntipnumber,"_",genetreenum,".pdf",sep="")
    pdf(file=file, width=5, height=5)
    
    par(mar=c(5.1, 4.1, 4.1, 4.1))
    plot(epsilonsdrealnew,border=NA,col=hcl.colors(10, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))
    # title("sd epsilons new")
    
    par(mar=c(5.1, 5.1, 5.1, 5.1))
    plot(heatmapdatanew,border=NA,col=hcl.colors(10, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))
    #title("average epsilons new")
    
    par(mar=c(5.1, 5.1, 5.1, 5.1))
    plot(spectreematrix,border=NA,col=hcl.colors(50, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))
    
    invisible(dev.off());
    
    save.image()
    
    setwd("..")
  }# close real mus loop
}#close real lambdasloop