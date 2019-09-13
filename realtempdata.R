t <- 0:100  # time
sig2 <- 0.01
## first, simulate a set of random deviates
x <- rnorm(n = length(t) - 1, sd = sqrt(sig2))
## now compute their cumulative sum
x <- c(30, cumsum(x)+30)
plot(t, x, type = "l", ylim = c(-2+30, 2+30))

real_temp_data=data.frame(t,x)
write.csv(real_temp_data,file="real_temp_data.csv")