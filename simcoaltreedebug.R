require("castor")
require("phybase")
require("ggtree")
require("treeAGG")
require("ape")


`simcoal` = function(rootnode,nodematrix,nspecies,seq,name)
  {
    theta<-nodematrix[rootnode,5]
    
    if(rootnode<=nspecies){
      {if(seq[rootnode] == 1){
        z<-list(gt="", height=as.matrix, node=as.matrix)
        z$gt<-name[rootnode]
        z$height<-0
        z$node<-nodematrix
        return(z)}
        else{
          treestr<-paste(name[rootnode],"s",1:seq[rootnode],sep="")
          i<-seq[rootnode]
          height<-rexp(1,rate=i*(i-1)/theta)
          brlens<-rep(0,i)
          father<-nodematrix[rootnode,1]
          fatherheight <- node.height(father,nodematrix,nspecies)
          
          while(height<fatherheight){
            nodematrix[rootnode,6]<-nodematrix[rootnode,6]+1
            ##randomly choose two nodes
            b<-sample(1:i,2)
            
            ##update groups
            newname<-paste("(",treestr[b[1]],sep="")
            newname<-paste(newname,":",sep="")
            newname<-paste(newname,round(height-brlens[b[1]],6),sep="")
            newname<-paste(newname,",",sep="")   
            newname<-paste(newname,treestr[b[2]],sep="")
            newname<-paste(newname,":",sep="")
            newname<-paste(newname,round(height-brlens[b[2]],6),sep="")
            newname<-paste(newname,")",sep="")
            treestr[b[1]]<-newname
            brlens[b[1]]<-height
            
            ##update dist,treestr, and branch length
            index<-1:i
            index[b[2]]<-0
            index<-index[index>0]
            treestr<-treestr[index]
            brlens<-brlens[index]
            if(i==2)
              break
            i<-i-1
            height<-height+rexp(1,rate=i*(i-1)/theta)
          }
          z<-list(gt="", height=as.matrix,node=as.matrix)
          z$gt<-treestr
          z$height<-brlens
          z$node<-nodematrix
          return(z)
        }
      }
      
    }
    if(rootnode>nspecies){
      son1<-nodematrix[rootnode,2]
      son2<-nodematrix[rootnode,3]
      leftstr<-simcoal(rootnode=son1,nodematrix=nodematrix,nspecies,seq,name)
      nodematrix<-leftstr$node
      rightstr<-simcoal(rootnode=son2,nodematrix=nodematrix,nspecies,seq,name)
      nodematrix<-rightstr$node
      i<-length(leftstr$gt)+length(rightstr$gt)
      
      treestr<-1:i
      treestr[1:length(leftstr$gt)]<-leftstr$gt
      treestr[(length(leftstr$gt)+1):i]<-rightstr$gt
      
      brlens<-1:i
      brlens[1:length(leftstr$height)]<-leftstr$height
      brlens[(length(leftstr$height)+1):i]<-rightstr$height
      
      
      
      height<-rexp(1,rate=i*(i-1)/theta)+ node.height(rootnode,nodematrix,nspecies)
      father<-nodematrix[rootnode,1]
      #edited below
      if(father == -9 | father == -8)
        fatherheight<-100000000000000000000000000000000000000 else
          fatherheight <- node.height(father,nodematrix,nspecies)
      brlens
      
      
      while(height<fatherheight){
        nodematrix[rootnode,6]<-nodematrix[rootnode,6]+1
        ##randomly choose two nodes
        b<-sample(1:i,2)
        
        ##update groups
        newname<-paste("(",treestr[b[1]],sep="")
        newname<-paste(newname,":",sep="")
        newname<-paste(newname,round(height-brlens[b[1]],6),sep="")
        newname<-paste(newname,",",sep="")   
        newname<-paste(newname,treestr[b[2]],sep="")
        newname<-paste(newname,":",sep="")
        newname<-paste(newname,round(height-brlens[b[2]],6),sep="")
        newname<-paste(newname,")",sep="")
        treestr[b[1]]<-newname
        brlens[b[1]]<-height
        
        ##update dist,treestr, and branch length
        index<-1:i
        index[b[2]]<-0
        index<-index[index>0]
        treestr<-treestr[index]
        brlens<-brlens[index]
        if(i==2)
          break
        i<-i-1
        height<-height+rexp(1,rate=i*(i-1)/theta)
      }
      if(nodematrix[rootnode,1]==-9 | nodematrix[rootnode,1]==-8)
        treestr<-paste(treestr,";",sep="")
      z<-list(gt="", height=as.matrix,node=as.matrix)
      z$gt<-treestr
      z$node<-nodematrix
      z$height<-brlens
      return(z)
    }
  }

## Not run:
# Generate a random tree with exponentially varying lambda & mu
# for(Ntips in c(rep(20,100),rep(100,100),rep(1000,100))){
Ntips=20
  rho = 1 # sampling fraction
  time_grid = seq(from=0, to=100, by=0.01)
  # for(lambdas in list(20+(100/tail(exp(0.1*time_grid),1))*exp(0.1*time_grid),0.2*time_grid+0.5,rep(2,length(time_grid)))){
  lambdas=0.2*time_grid+0.5
    mus = 0.1*time_grid
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
    
    # results = castor::get_clade_list(bigtree, postorder=TRUE, missing_value=-9)
    # nodematrix = list(nodes=cbind(results$clades, results$lengths,
    #                               matrix(-9,nrow=nrow(results$clades),ncol=3)),names=bigtree$tip.label, root=TRUE)
    # 
    genetreestuff = generate_gene_tree_msc(bigtree,allele_counts = 1,
                           population_sizes = runif(nspecies+Nnodes,min = 10^8,max=10^9),
                           generation_times = runif(nspecies+Nnodes,min = 0.001,max=1),
                           ploidy = 1);
    
    thetree=genetreestuff$tree
    
    ape::plot.phylo(thetree)
    title("genetree")
  
    # 
    # for(i in 1:nrow(nodematrix[["nodes"]])){
    #   nodematrix[["nodes"]][i,5]=runif(1,min=10^8,max=10^9)*(10^-11)*2
    #   #need to multiply above by mu (mutation rate per site
    #   #per generation) to get true theta
    # }
    # 
    # nspecies=length(bigtree[["tip.label"]])
    # 
    # rootnode=nrow(nodematrix[["nodes"]])
    # 
    # genetreestuff=simcoal(rootnode,nodematrix[["nodes"]],nspecies,seq=rep(1,nspecies),name=bigtree[["tip.label"]])
    # 
    # genetreeheight=genetreestuff[["height"]]
    # 
    # realgenetree=genetreestuff[["node"]]
    # 
    # genetreegt=genetreestuff[["gt"]]
    # 
    # thetree=castor::read_tree(genetreegt)
    # 
    # ape::plot.phylo(thetree)
    # title("genetree")
    
    # 
    # calculate true PDR, but this doens't make sense?
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
#   }#close functions loop
#   rm(list = ls())
# }#close Ntips loop
    