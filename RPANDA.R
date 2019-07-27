require("RPANDA")
require("diversitree")
require("BAMMtools")
require("castor")
require("mvMORPH")

real_temp_data=read.csv("real_temp_data.csv")

AICvals=vector()
t <- 0:100  # time
sig2 <- 0.01
## first, simulate a set of random deviates
x <- rnorm(n = length(t) - 1, sd = sqrt(sig2))
## now compute their cumulative sum
x <- c(30, cumsum(x)+30)
plot(t, x, type = "l", ylim = c(-2+30, 2+30))

random_data=data.frame(t,x)

f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.2, 0.01)
mu_par<-c()
result_exp <- sim_env_bd(real_temp_data,f.lamb,f.mu,lamb_par,mu_par,time.stop=10)
plot.phylo(result_exp[[1]])

real_fit=fit_env(phylo=result_exp[[1]],tot_time=max(node.age(result_exp[[1]])$ages),env_data=real_temp_data,f.lamb=f.lamb,f.mu=f.mu,lamb_par=lamb_par,mu_par=mu_par)

append(AICvals,as.numeric(real_fit['aicc']))

the_fit=fit_env(phylo=result_exp[[1]],tot_time=max(node.age(result_exp[[1]])$ages),env_data=random_data,f.lamb=f.lamb,f.mu=f.mu,lamb_par=lamb_par,mu_par=mu_par)

append(AICvals,as.numeric(the_fit['aicc']))

f.lamb <-function(t,x,y){y[1] + y[2] * x}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.2, 0.01)
mu_par<-c()

notreal_lin_fit=fit_env(phylo=result_exp[[1]],tot_time=max(node.age(result_exp[[1]])$ages),env_data=real_temp_data,f.lamb=f.lamb,f.mu=f.mu,lamb_par=lamb_par,mu_par=mu_par)

append(AICvals,as.numeric(notreal_lin_fit['aicc']))

# Fit the pure birth model (no extinction) with a constant speciation rate
f.lamb1 <-function(t,y){y[1]}
f.mu1<-function(t,y){0}
lamb_par1<-c(0.09)
mu_par1<-c()
result_cst <- fit_bd(phylo=result_exp[[1]],tot_time=max(node.age(result_exp[[1]])$ages),f.lamb1,f.mu1,lamb_par1,mu_par1,f=87/89,cst.lamb=TRUE,fix.mu=TRUE,dt=1e-3)

append(AICvals,as.numeric(result_cst['aicc']))

#do likelihood ratio test, but should I make one model 
#include y1,y2,y3,y4 for mu such that mu=y1exp(y2x)+y3+y4x? That way models are nested?

likerattest=-2*log(as.numeric(the_fit['LH'])-as.numeric(result_cst['LH']))

likerattest2=-2*log(as.numeric(real_fit['LH'])-as.numeric(result_cst['LH']))

likerattest3=-2*log(as.numeric(real_fit['LH'])-as.numeric(the_fit['LH']))

likerattest3=-2*log(as.numeric(real_fit['LH'])-as.numeric(the_fit['LH']))

likerattest4=-2*log(as.numeric(notreal_lin_fit['LH'])-as.numeric(real_fit['LH']))
#encapsulate all code in a loop that goes over increasing time.stop values in the simulated tree
#is this equivalent to increasing the amount of data in the tree? So LRT should be mor accurate?