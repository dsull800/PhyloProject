require("RPANDA")
require("diversitree")
require("castor")
require("phybase")
require("ggtree")
require("treeAGG")
require("ape")

t <- 0:10  # time
real_temp_data=data.frame(t,t)
f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.2, 0.01)
mu_par<-c()
result_exp <- sim_env_bd(real_temp_data,f.lamb,f.mu,lamb_par,mu_par,time.stop=100)
plot.phylo(result_exp[[1]])


#this code below gets the rootnode of the phylogeny but not actually
#mat=matTree(tree=result_exp[[1]
 
# sim.coaltree.sp(rootnode=mat[length(mat)],nodematrix=mat,nspecies=result_exp[["nblineages"]][length(result_exp[["nblineages"]])],seq=rep(1,result_exp[["nblineages"]][length(result_exp[["nblineages"]])]))

write.tree(phy=result_exp[[1]],file="newick1.txt")
# realtree=read.tree(file = "newick1")

#Be careful, the newickformat string ends in \n
fileName <- 'newick1.txt'
newickformat=readChar(fileName, file.info(fileName)$size)

nodematrix=read.tree.nodes(strsplit(newickformat,"\n"))

#I can set the population sizes by changing the entries of the 5th value of the matrix, 
#use a gamma distributed value or something, maybe.
for(i in 1:nrow(nodematrix$nodes)){
  nodematrix$nodes[i,5]=rgamma(n=1,shape=2,rate=1)
}

nspecies=tail(result_exp[[3]],n=1)

specnames=species.name(strsplit(newickformat,"\n"))


#nodematrix$nodes[nrow(nodematrix$nodes)-1,1]
rootnode=nrow(nodematrix$nodes)

genetreestuff=sim.coaltree.sp(rootnode,nodematrix$nodes,nspecies,seq=rep(1,nspecies),name=specnames)

genetreeheight=genetreestuff$height

realgenetree=genetreestuff$node

genetreegt=genetreestuff$gt

# plot.phylo(genetreegt)

thetree=read.tree(text=genetreegt)

plot.phylo(thetree)