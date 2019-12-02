require("castor")
require("prospectr")
require("plot.matrix")
require("matrixStats")
require("naniar")

funcnumber=0

wrkdir="/Users/danielsullivan/desktop/phylolambda0/PhyloProject-masterv"

# define plot colors
sim_color	= "blue"
emp_color 	= "red"
#set some loop variables for matrices and variables for grids
oldest_age_sim=1000
lineagecount=100
#set fineness of age grid
age_grid_fineness=.1
# specify number of extant tips in species tree
numtips=100000
#more varibales for column names
ncols=lineagecount
colnamesstuff=c()
for(i in seq(1,ncols)){
  inter_val=toString(i/ncols)
  colnamesstuff=c(colnamesstuff,inter_val)
}
# vector of R values
Rvec=c(10^-1,1,10)
#vector of rho values
rhovec=c(1,.5)
#number of spectrees to generate
numberofspec=7
#number of gene trees to generate for each species tree
numberofgen=10
#vector of lambda_0s
lambda_0vec=c(1,3,5)


#loop through functions for different scenarios
for(age2lambda in c(function(ages,lambda_0) rep(lambda_0,length(ages)),
                    function(ages,lambda_0) 0.1+ lambda_0*exp(-(ages-1)^2/(2*0.5^2)),
                    function(ages,lambda_0) 0.1 + lambda_0*exp(-0.005*ages),
                    function(ages,lambda_0) 0.7 + lambda_0*exp(0.003*ages))){
  
  
  for(age2mu in c(function(ages,lambda_0) rep(0,length(ages)), 
                  function(ages,lambda_0) 0.1 + 1*exp(-(ages-1)^2/(2*0.5^2)),
                  function(ages,lambda_0) 0.1 + 0.7*lambda_0*exp(-0.005*ages),
                  function(ages,lambda_0) 0.1 + 0.4*lambda_0*exp(0.0003*ages))){
    
    funcnumber=funcnumber+1
    
    for(rho in rhovec){
      
      spectreematrix=matrix(0,nrow=numberofspec,ncol=ncols)
      
      heatmapdatanew=matrix(0,nrow=length(lambda_0vec),ncol=ncols)
      rownames(heatmapdatanew)=lambda_0vec
      colnames(heatmapdatanew)=colnamesstuff
      
      heatmapdatanewsds=matrix(0,nrow=length(lambda_0vec),ncol=ncols)
      rownames(heatmapdatanew)=lambda_0vec
      colnames(heatmapdatanew)=colnamesstuff
      
      for(R in Rvec){
        
        #lambdasmatrix
        lambdamatrix=matrix(0,ncol=lineagecount,nrow=numberofgen*numberofspec)
        
        if(!dir.exists(paste(wrkdir,toString(funcnumber),"_",rho,"_",R,sep=""))){
          dir.create(paste(wrkdir,toString(funcnumber),"_",rho,"_",R,sep=""))
          setwd(paste(wrkdir,toString(funcnumber),"_",rho,"_",R,sep=""))
          dir.create("spectrees")
          dir.create("lambdaplots")
          dir.create("storedplots")
          dir.create("gentrees")
          dir.create("spectreeinfo")
          dir.create("matrixplots")
        }
        
        setwd(paste(wrkdir,toString(funcnumber),"_",rho,"_",R,sep=""))
        
        for(lambda_0 in lambda_0vec){
          
          #set working directory depending on loop variable
          
        #initialize loop variable
        Ntipnumber=-1
        for(Ntips in c(rep(numtips,numberofspec))){
          Ntipnumber=Ntipnumber+1
          
          #set age grid and rho
          age_grid_sim = seq(from=0, to = oldest_age_sim, by=age_grid_fineness)
          
          #simulate trees
          for(lambdas in list(age2lambda(age_grid_sim,lambda_0))){
            for(mus in list(age2mu(age_grid_sim,lambda_0))){
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
              cat(sprintf("Tree has %d tips, spans %g Myr, run number %g \n",length(spectree[["tip.label"]]),root_age,Ntipnumber))
              
              
              file=paste(lambda_0,"_",R,"_",Ntipnumber,"_",".txt",sep="")
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
              
              file=paste(Ntipnumber,"_",lambda_0,".pdf",sep="")
              pdf(file=file, width=5, height=5)
              plot(y=real_lambda_hat,x=lineagecountgrid,ylim=c(min(real_lambda_hat),max(real_lambda_hat)),type="l",col=emp_color)
              lines(y=lambda_hat_spectree,x=lineagecountgrid,type="l",col=sim_color)
              invisible(dev.off());
              
              file=paste(Ntipnumber,"_",lambda_0,".rds",sep="")
              saveRDS(c(real_lambda_hat,lambda_hat_spectree),file)
              
              setwd("..")
              
              for(i in 1:length(epsilonspectree)){
                spectreematrix[Ntipnumber,i]=spectreematrix[Ntipnumber,i]+epsilonspectree[i]
              }
              
              #want to simulate multiple genetrees for a given species tree
              genetreenum=-1
              
              gentime=10^-7
              popsize=spectreepdrpsr$PSR[1]/(R*gentime)
              for(i in seq(1,numberofgen)){
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
                
                #get lambdahatpprime values for later
                
                lambda_hat_p_prime_new=rev(gene_PSR)
                
                #compute epsilon values 
                
                epsilonnew=(lambda_hat_p_prime_new-real_lambda_hat)/real_lambda_hat
                
                # lambdamatrix stuff
                
                for(i in seq(1,length(epsilonnew))){
                  lambdamatrix[(Ntipnumber+1)*numberofgen-numberofgen+1+genetreenum,i]=epsilonnew[i]
                }
                
              }
                
                file=paste(Ntipnumber,"_",genetreenum,"_",R,"_",lambda_0,".pdf",sep="")
                
                # save plots of epsilon vs age and gene PSR
                setwd("storedplots")
                pdf(file=file, width=5, height=5)
                
                
                plot(y=epsilonnew,x=lineagecountgrid,type="l",col="green")
                title("realepsilon vs. time")
                
                plot(y=lambda_hat_p_prime_new,x=lineagecountgrid,type="l",col=emp_color)
                lines(y=real_lambda_hat,x=lineagecountgrid,type="l",col=sim_color)
                title("Gene PSR/Species PSR")
                legend(x=0, y=max(gene_PSR), legend=c("tree PSR", "calculated PSR"), lty=c(1,1), col=c(emp_color, sim_color))
                
                #this plot goes from past to present, reverse of what is standard in the rest of the code
                plot(gene_LTT$times-distancebetween, gene_LTT$lineages, type="l", xlab="time", ylab="# clades",col="red",ylim =c(0,101000))
                lines(lttcountspec$times, lttcountspec$lineages, type="l", xlab="time", ylab="# clades",ylim=c(0,101000))
                title("species tree/gene tree LTT")
                invisible(dev.off());
                
                # save information to files in certain directories
                file=paste(Ntipnumber,"_",genetreenum,"_",R,".txt",sep="")
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
        #close Ntips loop
        whatislambda_0=lambda_0vec==lambda_0
        for(i in seq(1:lineagecount)){
          heatmapdatanew[which(whatislambda_0),i]=colMeans(lambdamatrix)[i]
          
          heatmapdatanewsds[which(whatislambda_0),i]=colSds(lambdamatrix)[i]
        }
        }#close lambda_0 loop
      
      #create and plot matrices
      setwd("matrixplots")
      file=paste(Ntipnumber,"_",genetreenum,".pdf",sep="")
      pdf(file=file, width=5, height=5)
      
      par(mar=c(5.1, 4.1, 4.1, 4.1))
      plot(heatmapdatanewsds,border=NA,col=hcl.colors(10, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))
      # title("sd epsilons new")
      
      par(mar=c(5.1, 5.1, 5.1, 5.1))
      plot(heatmapdatanew,border=NA,col=hcl.colors(10, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))
      #title("average epsilons new")
      
      par(mar=c(5.1, 5.1, 5.1, 5.1))
      plot(spectreematrix,border=NA,col=hcl.colors(50, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE))
      
      invisible(dev.off());
      
      save.image()
      
      setwd("..")
        #close lambda_0 loop
      }#close R loop
      }#close rho loop
  }# close real mus loop
}#close real lambdasloop