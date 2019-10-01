X <- NULL #individual markov chain
X_n <- NULL #array of final chain values
X[1] <- 0 #markov chain values
X_n[1] <- 0
p <- 0.8 #probability
d <- NULL #random draw
n_chain_length <- 1000 #iterations
n_chains <- 100 #number of chains

g <- NULL #chain final value counter
i <- NULL #chain iteration counter

for (g in 1:(n_chains)) {

  for (i in 2:(n_chain_length)) {
    d <- runif(1)
    if (d <= p)  X[i] = X[i-1]+1
    else  X[i] = X[i-1]-1
  }
  
  X_n[g] <- X[n_chain_length]
}

results.vec <- NULL

results.vec['N'] <- length(X_n)
results.vec['mean'] <- mean(X_n)
results.vec['sd'] <- sd(X_n)
results.vec['median'] <- median(X_n)
results.vec['min'] <- min(X_n)
results.vec['max'] <- max(X_n)

hist(X_n)
results_p0.8 <- round(results.vec,1)