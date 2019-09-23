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













"simcoaltreespmu_mod"<-
  function(sptree, spname,seq,numgenetree,method="dirichlet",alpha=5.0)
  {
    nodematrix<-read.tree.nodes(sptree,spname)$nodes
    rootnode<-dim(nodematrix)[1]
    nspecies<-(rootnode+1)/2
    ntaxa<-sum(seq)
    
    for(i in 1:nrow(nodematrix)){
      nodematrix[i,5]=0.01
      #runif(1,min=10^8,max=10^9)
    }
    
    #generate mutation rates
    #nodematrix<-cbind(nodematrix,rep(-100,rootnode))
    if(tolower(method) == "gamma")
      nodematrix <- .mutation_exp(nodematrix,rootnode,rootnode,nspecies,alpha)
    if(tolower(method) == "dirichlet")
      nodematrix[,6] <- .rdirichlet(1,rep(alpha,dim(nodematrix)[1]))*dim(nodematrix)[1]
    if(tolower(method) == "user")
      nodematrix[,6] <- alpha
    
    index<-1
    seqname<-rep("",ntaxa)
    for(i in 1:nspecies)
      for(j in 1:seq[i])
      {
        if(seq[i] > 1)
          seqname[index]<-paste(spname[i],"s",j,sep="")
        else
          seqname[index]<-spname[i]
        index<-index+1
      }
    
    speciesmatrix<-matrix(0,nrow=nspecies,ncol=ntaxa)
    
    index<-1	
    for(i in 1:length(seq))
    {
      for(j in 1:seq[i])
      {
        speciesmatrix[i,index]<-1
        index<-index+1
      }
    }
    
    spnodedepth<-rep(0,2*nspecies-1)
    for(i in 1:(2*nspecies-1))
    {
      spnodedepth[i]<-node.height(i,nodematrix,nspecies)
    }
    
    treestr<-rep("",numgenetree)
    for(j in 1:numgenetree)
    {
      str<-sim.coaltree.sp(rootnode,nodematrix,nspecies,seq,name=spname)$gt
      genetree<-read.tree.nodes(str,name=seqname)$nodes
      genenodedepth<-rep(0,2*ntaxa-1)
      
      for(i in 1:(2*ntaxa-1))
      {
        genenodedepth[i]<-node.height(i,genetree,ntaxa)
      }	
      coaltree <- .populationMutation(nodematrix,spnodedepth,genetree,genenodedepth,speciesmatrix)
      treestr[j]<-write.subtree(dim(coaltree)[1],coaltree,seqname,dim(coaltree)[1])
    }
    z <- list(gt=as.character, st=as.matrix,seqname=as.character)
    z$gt <- treestr
    z$st <- nodematrix
    z$seqname<-seqname
    return(z)
  }

.rdirichlet <- function(n,a)
  ## pick n random deviates from the Dirichlet function with shape
  ## parameters a
{
  l<-length(a);
  x<-matrix(rgamma(l*n,a),ncol=l,byrow=TRUE);
  sm<-x%*%rep(1,l);
  x/as.vector(sm);
}

.mutation_exp<-function(sptree,root,inode,nspecies,alpha)
{
  if(inode == root)
  {
    sptree[root,6]<-1.0
    sptree <- .mutation_exp(sptree,root,sptree[root,2],nspecies,alpha)
    sptree <- .mutation_exp(sptree,root,sptree[root,3],nspecies,alpha)
  }
  else
  {
    sptree[inode,6]<-rgamma(1,alpha,(alpha/sptree[sptree[inode,1],6]))
    if(inode > nspecies)
    {
      sptree <- .mutation_exp(sptree,root,sptree[inode,2],nspecies,alpha)
      sptree <- .mutation_exp(sptree,root,sptree[inode,3],nspecies,alpha)
    }
  }
  return(sptree)
}


.populationMutation<-function (sptree, spnodedepth, genetree, genenodedepth, speciesmatrix)
{  
  index<-1
  
  nspecies<-(dim(sptree)[1]+1)/2
  ntaxa<-(dim(genetree)[1]+1)/2
  genetreenodes <- rep(-1,2*ntaxa)
  
  for(i in 1:nspecies)
  {
    seq<-(1:ntaxa)*(speciesmatrix[i,])
    seq<-seq[seq>0]
    for(inode in 1:length(seq))
    {        
      inodegene <- seq[inode];
      stop=0;
      while(inodegene != dim(genetree)[1])
      {
        #check if the node is already taken care of
        for(k in 1:index)
          if(inodegene == genetreenodes[k]) 
          {
            stop<-1
            break
          }
        if(stop == 1) 	break
        
        #change the branch length of node p		
        genetree[inodegene,4] <- .ChangeBrlen(sptree, spnodedepth, i, genetree, genenodedepth, inodegene)
        
        #copy p to genetreenode
        genetreenodes[index] <- inodegene
        index<-index+1
        
        #reset inodegene
        inodegene<-genetree[inodegene,1]
      }
    }
  }
  genetree[,4]<-round(genetree[,4],6)
  return (genetree)
  
}

.ChangeBrlen<-function(sptree, spnodedepth, spnode, genetree, genenodedepth, genenode)
{
  inode <- .FindSpnodeDownGenenode(sptree,spnodedepth,spnode,genenodedepth,genenode)
  jnode <- .FindSpnodeDownGenenode(sptree,spnodedepth,spnode,genenodedepth,genetree[genenode,1])
  
  if(inode == jnode)
  {
    length <- (genetree[genenode,4]) * sptree[inode,6]
  }
  else
  {
    father <- sptree[inode,1]
    length <- (spnodedepth[father] - genenodedepth[genenode])*sptree[inode,6]
    while(father != jnode)
    {
      inode <- father;
      father <- sptree[father,1] 
      length <- length + (spnodedepth[father] - spnodedepth[inode])*(sptree[inode,6])
    }
    length <- length + (genenodedepth[genetree[genenode,1]] - spnodedepth[father])*sptree[father,6]
  }
  return(length)		
}

.FindSpnodeDownGenenode<-function(sptree, spnodedepth, spnode, genenodedepth, genenode)
{
  findnode<-spnode;
  root<-dim(sptree)[1]
  
  depth <- genenodedepth[genenode]
  father <- sptree[spnode,1]
  
  if(genenode > (length(genenodedepth)+1)/2)
  {
    while(spnodedepth[father] <= depth)
    {
      if(father == root)
      {
        findnode <- father
        break
      }
      else
      {
        findnode <- father
        father <- sptree[father,1]
      }
    }
  }    
  return (findnode)
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

# results = castor::get_clade_list(bigtree, postorder=TRUE, missing_value=-9)
# nodematrix = list(nodes=cbind(results$clades, results$lengths,
#                               matrix(-9,nrow=nrow(results$clades),ncol=3)),names=bigtree$tip.label, root=TRUE)

# 
# for(i in 1:nrow(nodematrix[["nodes"]])){
#   nodematrix[["nodes"]][i,5]=runif(1,min=10^8,max=10^9)*(10^-11)*2
#   #need to multiply above by mu (mutation rate per site
#   #per generation) to get true theta
# }

nspecies=length(bigtree[["tip.label"]])

# rootnode=nrow(nodematrix[["nodes"]])
# 
bigtreestring=write_tree(bigtree)

noclock=simcoaltreespmu_mod(bigtreestring,spname=bigtree[["tip.label"]],numgenetree=1,alpha=5,method="dirichlet",seq=rep(1,nspecies))


genetree=read_tree(noclock$gt)
plot.phylo(genetree)
title("genetree")

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
# # 
# # calculate true PDR, but this doens't make sense?
# lambda_slopes = diff(lambdas)/diff(time_grid);
# lambda_slopes = c(lambda_slopes[1],lambda_slopes)
# PDRs = lambdas - mus - (lambda_slopes/lambdas)
# # Fit PDR on grid
# Ngrid = 10
# height=max(get_all_distances_to_root(thetree))
# age_grid = seq(0,height,length.out=Ngrid)
# # ERROR: Fitting failed: Provided age-grid range (0 - 1.20583) does not cover entire required age range (0 - 1.20583) look at lines 50-51 in fit_hbd_pdr_on_grid
# fitpdr = castor::fit_hbd_pdr_on_grid(thetree,
#                                      age_grid=age_grid,
#                                      min_PDR = -20,
#                                      max_PDR = +50,
#                                      condition = "crown",
#                                      Ntrials = 10,# perform 10 fitting trials
#                                      Nthreads = 4,# use two CPUs
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
#         xlim = c(genetreeheight,0),
#         ylim=c(-100,100))
#   # xlim = c(-100,100))
#   lines(x = seq(from=genetreeheight,to=0,length.out=length(time_grid)),
#         y = PDRs,
#         type = 'l',
#         col = 'blue');
# }
# 
# # Fit PSR on grid
# oldest_age=genetreeheight/2 # only consider recent times when fitting
# Ngrid = Ntips
# age_grid = seq(from=0,to=oldest_age,length.out=Ngrid)
# fit = fit_hbd_psr_on_grid(thetree,
#                           oldest_age = oldest_age,
#                           age_grid = age_grid,
#                           min_PSR = -50,
#                           max_PSR = +50,
#                           guess_PSR=1,
#                           # fixed_PSR=rep(1,Ngrid),
#                           condition = "crown",
#                           Ntrials = 10,# perform 10 fitting trials
#                           Nthreads = 5,# use two CPUs
#                           max_model_runtime = 1) # limit model evaluation to 1 second
# if(!fit[["success"]]){
#   cat(sprintf("ERROR: Fitting failed: %s\n",fit[["error"]]))
# }else{
#   cat(sprintf("Fitting succeeded:\nLoglikelihood=%g\n",fit[["loglikelihood"]]))
#   # plot fitted PSR
#   plot( x = fit[["age_grid"]],
#         y = fit[["fitted_PSR"]],
#         main = 'Fitted PSR',
#         xlab = 'age',
#         ylab = 'PSR',
#         type = 'b',
#         xlim = c(genetreeheight/2,0))
#   # plot deterministic LTT of fitted model
#   plot( x = fit[["age_grid"]],
#         y = fit[["fitted_LTT"]],
#         main = 'Fitted dLTT',
#         xlab = 'age',
#         ylab = 'lineages',
#         type = 'b',
#         log = 'y',
#         xlim = c(genetreeheight/2,0))
# }
# lttcount=castor::count_lineages_through_time(thetree,Ntimes=100)
# plot(lttcount$times, lttcount$lineages, type="l", xlab="time", ylab="# clades")
#   }#close functions loop
#   rm(list = ls())