require("RPANDA")
require("diversitree")
require("BAMMtools")
require("castor")
require("mvMORPH")

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
result_exp <- sim_env_bd(random_data,f.lamb,f.mu,lamb_par,mu_par,time.stop=10)
plot.phylo(result_exp[[1]])
the_fit=fit_env(phylo=result_exp[[1]],tot_time=max(node.age(result_exp[[1]])$ages),env_data=random_data,f.lamb=f.lamb,f.mu=f.mu,lamb_par=lamb_par,mu_par=mu_par)
