library(MASS) # Generation of samples from multivariate normal distributions
library(corpcor) # Ledoit Wolf shrinkage estimator
library(ddpcr) # Gets cov.shrink() to STFU
library(purrr) # calculating dot products by mapping sample matrices to covariance
library(ggplot2) # customised graphics
library(gridExtra) # arranging arrays of ggplots
library(reshape2)
library(stringr)

## Preliminary: generate samples from a multivariate normal distribution with mean vector = 0 in each element and specified covariance matrices
## Covariance matrices take 4 forms: identity, identity linear decaying, identity exponentially decaying and toeplitz
## For each covariance structure, we generate samples of a size dependant on the dimensionality (the length, or number of rows of the matrix)
## So, sample size = p, sample size = p*log(p), sample size = 5*p*log(p) and sample size = 10*p*log(p)
## P, the dimensionality (read: matrix width) of the samples: 5,10,50,100,250
## For each combination of covariance structure, sample size and dimensionality, we repeat the simulation 80 times, using replicate()
## This means we generate 5 (no. of sample sizes)*5 (no of dimensionalities)*4 (number of covariance matrices) * 20 (number of replications)
## Hence, we generate 2,000 matrices

#---- Part 1: Matrix Generation  ----

# experimental setup - dimensionality, mean vectors (all 0), sample sizes (a function of dimensionality)
p = c(5,10,50,100,250)
mean_vectors = sapply(X = p,FUN = rep,x = 0, times = p)
sample_size = list(p,ceiling(p*log(p)),ceiling(5*p*log(p)),ceiling(10*p*log(p)))

# identity matrix generation
# example
diag(3)

# linear decay matrix generation
diag_linear_decay <- function(n) {
	return(diag(n:1,n,n))
}

# example 
diag_linear_decay(3)

# exponential decay matrix generation
exp_decay <- function(n) {
		v = c()
		for (i in 1:n) {
			v[i] <- n^(2-i)
		}
		return(diag(v,n,n))
}

# example 
exp_decay(3)

# toeplitz matrix generation
toeplitz <- function(n) {
	i = 1
	j = 1
	m = matrix(nrow = n,ncol = n)
	for (j in 1:n) {
		for (i in 1:n) {
			m[i,j] = 0.5^(abs(i-j))
		}
	}
	return(m)
}

# example
toeplitz(3)

# matrix generation function
gen_samples <- function(n,mean_vectors,p,sigma_type) {
	
	res = switch(sigma_type,list("identity" = list()),list("linear_decay" = list()),list("exp_decay" = list()),list("toeplitz" = list()))
	
	for (i in 1:length((sample_size))) {
		if (sigma_type == 1) {
				res[[1]][[i]] = mapply(mvrnorm,n = sample_size[[i]], mu = mean_vectors, Sigma = sapply(p,diag))
			} else if (sigma_type == 2) {
				res[[1]][[i]] = mapply(mvrnorm,n = sample_size[[i]], mu = mean_vectors, Sigma = sapply(p,diag_linear_decay))
			} else if (sigma_type == 3) {
				res[[1]][[i]] = mapply(mvrnorm,n = sample_size[[i]], mu = mean_vectors, Sigma = sapply(p,exp_decay))	
			} else if (sigma_type == 4) {	
				res[[1]][[i]] = mapply(mvrnorm,n = sample_size[[i]], mu = mean_vectors, Sigma = sapply(p,toeplitz))	
			}
		names(res[[1]][[i]]) <- c("p=5","p=10","p=50","p=100","p=250")
	}
	names(res[[1]]) <- c("n=p","n=p*log(p)", "n=5p*log(p)","n=10p*log(p)")
	return(res)
}

# generation of matrices (samples from multivariate normal with specified covariance structure)
n_rep <- 20

res <- gen_samples(n = sample_size,mean_vectors = mean_vectors,p = p,sigma_type = 1)
res_rep_sigma1 <- replicate(n_rep, gen_samples(n = sample_size,mean_vectors = mean_vectors,p = p,sigma_type = 1))
res_rep_sigma2 <- replicate(n_rep, gen_samples(n = sample_size,mean_vectors = mean_vectors,p = p,sigma_type = 2))
res_rep_sigma3 <- replicate(n_rep, gen_samples(n = sample_size,mean_vectors = mean_vectors,p = p,sigma_type = 3))
res_rep_sigma4 <- replicate(n_rep, gen_samples(n = sample_size,mean_vectors = mean_vectors,p = p,sigma_type = 4))

#---- Part 2: Covariance estimation ----

## The objective now is to investigate the estimation of the population covariance structure of the independant variables from which these matrices are sampled
## The population covariance matrix structure is known
## Here, we use three types of estimators: maximum likelihood with known mean, a diagonal estimator of only variances, and the Ledoit-Wolf shrinkage estimator

ml_cov <- function(m) {
# classic maximum likelihood estimator
# as we know the mean vector = 0 in each element, this drops out from the calculation
(1/nrow(m)) * (t(m) %*% m)
}

# we see that the maximum likelihood estimate of covariance is very similar to the inbuilt R estimate (which would, in contrast, use an estimate of the mean vector)
ml_cov(res_rep_sigma1$identity$`n=10p*log(p)`$`p=5`)
cov(res_rep_sigma1$identity$`n=10p*log(p)`$`p=5`)

# leidot-wolf estimate
cov.shrink(res_rep_sigma1$identity$`n=p`$`p=5`)

# We check that the package is doing things correctly by translating the matlab code;
LedoitWolf_cov_estimate <- function(x,shrink,return_shrink = F) {
  # demean matrix
  size_x <- dim(x)
  t <- size_x[1] 
  n <- size_x[2] 
  meanx <- apply(FUN = mean,x,MARGIN = 2)
  x <- x-t(replicate(nrow(x),apply(FUN = mean,x,MARGIN = 2))) # as per matlab code
  
  # sample covariance matrix - assuming mean vector = 0 ???
  sample_cov <- (1/t)*(t(x) %*% x)
  #print(sample_cov) #debugging
  
  # compute prior
  mean_var <- mean(diag(sample_cov)) # mean of diagonal entries of sample covariance matrix
  prior <- mean_var * diag(1,n)
  	
  	# compute shrinkage parameters if none provided
  	if (missing(shrink)) {
  	# what they call p
  	y <- x^2
  	phiMat <- t(y) %*% y/t-sample_cov^2
  	phi <- sum(phiMat)
  	# what they call c
  	gamma <- norm(sample_cov - prior,type = "F")^2 # Frobenius Norm
  	
  	# compute shrinkage constant
  	kappa <- phi/gamma
  	shrink <- max(0,min(1,kappa/t))
  		
  	# checking atomicity and length is the standard scalar check in R, as far as I am aware, believe it or not
  	} else if (length(shrink) != 1 || !is.atomic(shrink)) { # nargin = number of function input arguments
  		stop("shrink must be atomic and of length 1")
  	}
  
  sigma <- shrink * prior+(1-shrink) * sample_cov
  if (!return_shrink) {
    return(sigma)
  } else {
    return(list('sigma' = sigma,'shrinkage' = shrink))
  	#beep() #  )': no ding
  }

}

# checking package functions against custom function translated and tested against matlab code
x <- t(matrix(1:9,nrow = 3))
linshrink_cov(x) # differs from matlab implementation
cov.shrink(x) # differs from matlab implementation
LedoitWolf_cov_estimate(x,return_shrink = T)

# only translated code works identically to matlab function (verified using home copy of matlab)

# Finally, we simply use the diagonal entries of the ML estimate as our 'diagonal estimator'
# So we simply create a wrapper function for the existing ML estimator
diag_ml_cov <- function(x) {
	return(diag(diag(ml_cov(x)))) # nested diag so it returns covariance as a matrix, not a vector
}

# example of each covariance estimator
print(x)
ml_cov(x)
diag_ml_cov(x)
LedoitWolf_cov_estimate(x,return_shrink = T)

matrices <- list(`identity` = res_rep_sigma1,
                 `linear` = res_rep_sigma2,
                 `expon` = res_rep_sigma3,
                 `toeplitz` = res_rep_sigma4)

# generating covariance matrice estimates for each matrix, using each type of estimator - this results in 6,000 matrices 
ml_cov_estimates <- rapply(matrices, ml_cov, how="replace")
diag_ml_cov_estimates <- rapply(matrices, diag_ml_cov, how="replace")
#LedoitWolf_cov_estimates <- rapply(matrices, LedoitWolf_cov_estimate, how="replace",return_shrink = F) 
# due to computational time involved, estimating on a seperate machine - FAILED
# works (slowly on individual matrices) but fails when recursively applied over list...
# due to computational time - had to use package version (cov.shrink)
# quiet(LedoitWolf_cov_estimates <- rapply(matrices, cov.shrink, how="replace")) # WAY faster
# LW Estimates calculated on laptop overnight - saved as RDS, imported here
LedoitWolf_cov_estimates <- readRDS("lw_est.rds")



#---- Part 3: Assessing Covariance Estimates  ----

gen_covariance_structures <- function(p,sigma_type) {
	
	res = switch(sigma_type,list("identity" = list()),list("linear_decay" = list()),list("exp_decay" = list()),list("toeplitz" = list()))
	
		if (sigma_type == 1) {
			res[[1]] = sapply(p,diag)
		} else if (sigma_type == 2) {
			res[[1]] = sapply(p,diag_linear_decay)
		} else if (sigma_type == 3) {
			res[[1]] = sapply(p,exp_decay)
		} else if (sigma_type == 4) {	
			res[[1]] = sapply(p,toeplitz)
		names(res[[1]]) <- c("p=5","p=10","p=50","p=100","p=250")
	}
	return(res)
}

identity_cov_structure <- gen_covariance_structures(p=p,sigma = 1)
linear_cov_structure <- gen_covariance_structures(p=p,sigma = 2)
expon_cov_structure <- gen_covariance_structures(p=p,sigma = 3)
toeplitz_cov_structure <- gen_covariance_structures(p=p,sigma = 4)

cov_structures <- list(`identity` = identity_cov_structure,
											 `linear` = linear_cov_structure,
											 `exponential` = expon_cov_structure,
											 `toeplitz` = toeplitz_cov_structure)

#---- Question 1: How well are the individual eigenvalues (variances) estimated? ----

# Estimation of Eigenvalues (variances)
# Largest and Smallest Eigenvalues compared with groundtruth 

eigen_range <- function(x) {
	# returns largest and smallest eigenvalues
	return(range(eigen(x,only.values = T,symmetric = T)))
}

# Largest / Smallest Estimates

eigen_range_ml <- rapply(ml_cov_estimates,eigen_range, how="replace")
eigen_range_diag_ml <- rapply(diag_ml_cov_estimates,eigen_range, how="replace")
eigen_range_LedoitWolf <- rapply(LedoitWolf_cov_estimates,eigen_range, how="replace")

# Largest / Smallest Ground Truth

eigen_range_cov_structures <- rapply(cov_structures,eigen_range, how="replace")

#---- Question 2: How well is the trace of the covariance matrix estimated? ----

# Estimation of trace of covariance matrix
# Sum of diagonal entries (since we know the matrix is always square)

tr <- function(x) { sum(diag(x)) }

# Trace of Estimates

trace_ml <- rapply(ml_cov_estimates,tr, how="replace")
trace_diag_ml <- rapply(diag_ml_cov_estimates,tr, how="replace")
trace_LedoitWolf <- rapply(LedoitWolf_cov_estimates,tr, how="replace")

# Trace of Ground Truth

trace_cov_structures <- rapply(cov_structures,tr, how="replace")



#---- Question 3: How well are the principal components (eigenvectors) of the covariance matrix estimated? ----

# Just need to look at first principal component (as if this is wrong, so are others)
# In R the output of eigen is normalised anyway - so the denominator of the normalised dot product will be one, and can hence be left out

eigen_dot_product <- function(x,m) {
		return(eigen(x)$vec[,1] %*% eigen(m)$vec[,1])
	}

# Normalised Dot Products of Estimates (x) vs Ground Truth (m)

	
# https://stackoverflow.com/questions/49252400/r-purrr-flatten-list-of-named-lists-to-list-and-keep-names
my_flatten <- function (x, use.names = TRUE, classes = "ANY") 
{
	#' Source taken from rlist::list.flatten
	len <- sum(rapply(x, function(x) 1L, classes = classes))
	y <- vector("list", len)
	i <- 0L
	items <- rapply(x, function(x) {
		i <<- i + 1L
		y[[i]] <<- x
		TRUE
	}, classes = classes)
	if (use.names && !is.null(nm <- names(items))) 
		names(y) <- nm
	y
}

 #mapply(FUN = eigen_dot_product,x = ml_cov_estimates[[i]],m = cov_structures)
 
# ML
 identity_cov_flat <- my_flatten(ml_cov_estimates[[1]])
 identity_cov_struct_flat <- rep(my_flatten(cov_structures[[1]]),80)
 dot_products_ml_identity <- map2(identity_cov_flat, identity_cov_struct_flat, eigen_dot_product)
 
 linear_cov_flat <- my_flatten(ml_cov_estimates[[2]])
 linear_cov_struct_flat <- rep(my_flatten(cov_structures[[2]]),80)
 dot_products_ml_linear <- map2(identity_cov_flat, identity_cov_struct_flat, eigen_dot_product)
 
 exponential_cov_flat <- my_flatten(ml_cov_estimates[[3]])
 exponential_cov_struct_flat <- rep(my_flatten(cov_structures[[3]]),80)
 dot_products_ml_exponential <- map2(exponential_cov_flat, exponential_cov_struct_flat, eigen_dot_product)
 
 toeplitz_cov_flat <- my_flatten(ml_cov_estimates[[4]])
 toeplitz_cov_struct_flat <- rep(my_flatten(cov_structures[[4]]),80)
 dot_products_ml_toeplitz <- map2(toeplitz_cov_flat, toeplitz_cov_struct_flat, eigen_dot_product)
 
 # Diag ML
 identity_cov_flat <- my_flatten(diag_ml_cov_estimates[[1]])
 identity_cov_struct_flat <- rep(my_flatten(cov_structures[[1]]),80)
 dot_products_diag_ml_identity <- map2(identity_cov_flat, identity_cov_struct_flat, eigen_dot_product)
 
 linear_cov_flat <- my_flatten(diag_ml_cov_estimates[[2]])
 linear_cov_struct_flat <- rep(my_flatten(cov_structures[[2]]),80)
 dot_products_diag_ml_linear <- map2(identity_cov_flat, identity_cov_struct_flat, eigen_dot_product)
 
 exponential_cov_flat <- my_flatten(diag_ml_cov_estimates[[3]])
 exponential_cov_struct_flat <- rep(my_flatten(cov_structures[[3]]),80)
 dot_products_diag_ml_exponential <- map2(exponential_cov_flat, exponential_cov_struct_flat, eigen_dot_product)
 
 toeplitz_cov_flat <- my_flatten(diag_ml_cov_estimates[[4]])
 toeplitz_cov_struct_flat <- rep(my_flatten(cov_structures[[4]]),80)
 dot_products_diag_ml_toeplitz <- map2(toeplitz_cov_flat, toeplitz_cov_struct_flat, eigen_dot_product)
 
 # Ledoit-Wolf
 identity_cov_flat <- my_flatten(LedoitWolf_cov_estimates[[1]])
 identity_cov_struct_flat <- rep(my_flatten(cov_structures[[1]]),80)
 dot_products_LedoitWolf_cov_estimates_identity <- map2(identity_cov_flat, identity_cov_struct_flat, eigen_dot_product)
 
 linear_cov_flat <- my_flatten(LedoitWolf_cov_estimates[[2]])
 linear_cov_struct_flat <- rep(my_flatten(cov_structures[[2]]),80)
 dot_products_LedoitWolf_cov_estimates_linear <- map2(identity_cov_flat, identity_cov_struct_flat, eigen_dot_product)
 
 exponential_cov_flat <- my_flatten(LedoitWolf_cov_estimates[[3]])
 exponential_cov_struct_flat <- rep(my_flatten(cov_structures[[3]]),80)
 dot_products_LedoitWolf_cov_estimates_exponential <- map2(exponential_cov_flat, exponential_cov_struct_flat, eigen_dot_product)
 
 toeplitz_cov_flat <- my_flatten(LedoitWolf_cov_estimates[[4]])
 toeplitz_cov_struct_flat <- rep(my_flatten(cov_structures[[4]]),80)
 dot_products_LedoitWolf_cov_estimates_toeplitz <- map2(toeplitz_cov_flat, toeplitz_cov_struct_flat, eigen_dot_product)
 
 #---- Question 4: How well is the covariance matrix estimated overall? ----

 # For example you could look at the sum of squared errors in Σ−Σˆ which is the trace of (Σ−Σˆ )(Σ−Σˆ ).
 
 SSE_cov <- function(x,m) {
 	return(tr((x-m) %*% (x-m)))
 }
 
 # ML
 identity_cov_flat <- my_flatten(ml_cov_estimates[[1]])
 identity_cov_struct_flat <- rep(my_flatten(cov_structures[[1]]),80)
 SSE_ml_identity <- map2(identity_cov_flat, identity_cov_struct_flat, SSE_cov)
 
 linear_cov_flat <- my_flatten(ml_cov_estimates[[2]])
 linear_cov_struct_flat <- rep(my_flatten(cov_structures[[2]]),80)
 SSE_ml_linear <- map2(identity_cov_flat, identity_cov_struct_flat, SSE_cov)
 
 exponential_cov_flat <- my_flatten(ml_cov_estimates[[3]])
 exponential_cov_struct_flat <- rep(my_flatten(cov_structures[[3]]),80)
 SSE_ml_exponential <- map2(exponential_cov_flat, exponential_cov_struct_flat, SSE_cov)
 
 toeplitz_cov_flat <- my_flatten(ml_cov_estimates[[4]])
 toeplitz_cov_struct_flat <- rep(my_flatten(cov_structures[[4]]),80)
 SSE_ml_toeplitz <- map2(toeplitz_cov_flat, toeplitz_cov_struct_flat, SSE_cov)
 
 # Diag ML
 identity_cov_flat <- my_flatten(diag_ml_cov_estimates[[1]])
 identity_cov_struct_flat <- rep(my_flatten(cov_structures[[1]]),80)
 SSE_diag_ml_identity <- map2(identity_cov_flat, identity_cov_struct_flat, SSE_cov)
 
 linear_cov_flat <- my_flatten(diag_ml_cov_estimates[[2]])
 linear_cov_struct_flat <- rep(my_flatten(cov_structures[[2]]),80)
 SSE_diag_ml_linear <- map2(identity_cov_flat, identity_cov_struct_flat, SSE_cov)
 
 exponential_cov_flat <- my_flatten(diag_ml_cov_estimates[[3]])
 exponential_cov_struct_flat <- rep(my_flatten(cov_structures[[3]]),80)
 SSE_diag_ml_exponential <- map2(exponential_cov_flat, exponential_cov_struct_flat, SSE_cov)
 
 toeplitz_cov_flat <- my_flatten(diag_ml_cov_estimates[[4]])
 toeplitz_cov_struct_flat <- rep(my_flatten(cov_structures[[4]]),80)
 SSE_diag_ml_toeplitz <- map2(toeplitz_cov_flat, toeplitz_cov_struct_flat, SSE_cov)
 
 # Ledoit-Wolf
 identity_cov_flat <- my_flatten(LedoitWolf_cov_estimates[[1]])
 identity_cov_struct_flat <- rep(my_flatten(cov_structures[[1]]),80)
 SSE_LedoitWolf_cov_estimates_identity <- map2(identity_cov_flat, identity_cov_struct_flat, SSE_cov)
 
 linear_cov_flat <- my_flatten(LedoitWolf_cov_estimates[[2]])
 linear_cov_struct_flat <- rep(my_flatten(cov_structures[[2]]),80)
 SSE_LedoitWolf_cov_estimates_linear <- map2(identity_cov_flat, identity_cov_struct_flat, SSE_cov)
 
 exponential_cov_flat <- my_flatten(LedoitWolf_cov_estimates[[3]])
 exponential_cov_struct_flat <- rep(my_flatten(cov_structures[[3]]),80)
 SSE_LedoitWolf_cov_estimates_exponential <- map2(exponential_cov_flat, exponential_cov_struct_flat, SSE_cov)
 
 toeplitz_cov_flat <- my_flatten(LedoitWolf_cov_estimates[[4]])
 toeplitz_cov_struct_flat <- rep(my_flatten(cov_structures[[4]]),80)
 SSE_LedoitWolf_cov_estimates_toeplitz <- map2(toeplitz_cov_flat, toeplitz_cov_struct_flat, SSE_cov)
 
 #---- Question 5: How stable are the estimates? ----
 
 
 
 #---- Plotting ----
 
 # ML,D-ML,LW Plotting Function
plot_eigens <- function(x,y,z,start,end,covariance_structure,density = F) {
	options(digits = 3, scipen = -2)
	# x is estimator list of eigens
	 eigen_smallest <- list()
	 eigen_biggest <- list()
	 df_small <- data.frame(matrix(NA, ncol=1, nrow=21)[-1]) # from SO
	 df_large <- data.frame(matrix(NA, ncol=1, nrow=21)[-1]) # from SO
	 for (i in 1:20) {
	 	title <- names(my_flatten(x)[start+i])
	 	title <- str_remove(pattern = str_extract(pattern = "[^.]*.[^.]*",string = title),string = title)
	 	#title <- substr(title,17,title) # trim off "identity.identity"
	 	df_x <- as.data.frame(do.call(rbind, my_flatten(x)[seq(start+i,end,by = 20)]))
	 	df_y <- as.data.frame(do.call(rbind, my_flatten(y)[seq(start+i,end,by = 20)]))
	 	df_z <- as.data.frame(do.call(rbind, my_flatten(z)[seq(start+i,end,by = 20)]))
	 	colnames(df_x) <- c('Smallest_ML',"Largest_ML")
	 	colnames(df_y) <- c('Smallest_D_ML',"Largest_D_ML")
	 	colnames(df_z) <- c('Smallest_TP',"Largest_TP")
	 	quiet(df <- melt(cbind(df_x,df_y,df_z)))
	 	df_small <- df[(df[,1] == "Smallest_ML" | df[,1] == "Smallest_D_ML" | df[,1] == "Smallest_TP" ),]
	 	df_large <- df[(df[,1] == "Largest_ML" | df[,1] == "Largest_D_ML" | df[,1] == "Largest_TP" ),]
		if (density == T) {
		
			g <- ggplot(data = df_small) + 
				geom_density(mapping = aes(x = value,color = variable,fill = variable),alpha = 0.5) +
				geom_vline(xintercept = 1,color = 'red') +
				labs(title = title,y = "density", x = "") #+
				theme_light() +
				theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 0.5,size = 6)) +
				theme(plot.title = element_text(size = 5)) + theme(legend.position="none") +
				scale_x_continuous(breaks = NULL)
			
			} else if (covariance_structure == "Identity") {
	 		g <- ggplot(data = df_small) + 
	 			geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
				geom_vline(xintercept = 1,color = 'red') +
	 			labs(title = title,y = "count", x = "") +
	 			theme_light() +
	 			theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 0.5,size = 6)) +
	 			theme(plot.title = element_text(size = 5)) + theme(legend.position="none")  #+	
	 	} else if (covariance_structure == "Exponential") {
		 	g <- ggplot(data = df_small) + 
		 		geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
		 		labs(title = title,y = "count", x = "") +
		 		theme_light() +
		 		theme(plot.title = element_text(size = 5)) +
		 		theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 0.5,size = 6)) + 
				theme(legend.position="none")
	 	} else {
	 		g <- ggplot(data = df_small) + 
	 			geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
	 			labs(title = title,y = "count", x = "") +
	 			theme_light() +
	 			theme(plot.title = element_text(size = 5)) + theme(legend.position="none") #+
	 	}
	 	
	 	if (covariance_structure == "Identity") {
	 		gl <- ggplot(data = df_large) + 
	 			geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
	 			geom_vline(xintercept = 1,color = 'red') +
	 			labs(title = title,y = "count", x = "") +
	 			theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 0.5,size = 6)) +
				theme_light() +
	 			theme(plot.title = element_text(size = 5)) + theme(legend.position="none") #+	
	 	} else if (covariance_structure == "Exponential") {
	 		gl <- ggplot(data = df_large) + 
	 			geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
	 			labs(title = title,y = "count", x = "") +
	 			theme_light() +
	 			theme(plot.title = element_text(size = 5)) +
	 			theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 0.5,size = 6)) + 
				theme(legend.position="none") 
	 	} else {
	 		gl <- ggplot(data = df_large) + 
	 			geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
	 			labs(title = title,y = "count", x = "") +
	 			theme_light() +
	 			theme(plot.title = element_text(size = 5)) + theme(legend.position="none") #+
	 	}
	 	
	 	#colnames(df) <- c(paste('Smallest_',title),paste("Largest_",title))
#
	 	#df_small[,i] <- df[,1]
	 	#df_large[,i] <- df[,2]
	 	#colnames(df_small)[i] <- title
	 	#colnames(df_large)[i] <- title
	 	
	 	eigen_smallest[[i]] <- g
	 	eigen_biggest[[i]] <- gl
	 }
	 density_ <- NULL;dotplot_ <- NULL
	 if (density) { density_ <- "Density"} else { dotplot_ <- "Dotplots" }
	 
	 quiet(grob <- grid.arrange(
	 	grobs = eigen_smallest,
	 	ncol = 5, 
	 	top = paste(covariance_structure," Population Covariance Structure"),
	 	bottom = paste("Estimators ", "Smallest Eigenvalues")
	 ))
	 quiet(grob_large <- grid.arrange(
	 	grobs = eigen_biggest,
	 	ncol = 5, 
	 	top = paste(covariance_structure," Population Covariance Structure"),
	 	bottom = paste("Estimators ", "Largest Eigenvalues - ",dotplot_,density_)
	 ))
	 return(list('small' = df_small,'large' = df_large,'g' = grob,'gl' = grob_large))
}
	# we expect largest and smallest values to be 1
for (i in 1:2) {
	if (i == 1) { 
		bool = F
		density = ""
		dotplot = "dotplots"
	}	else {
		 bool = T
		 density = "density"
		 dotplot = ""
	}
	 eig_ident.ls <- plot_eigens(eigen_range_ml,eigen_range_diag_ml,eigen_range_LedoitWolf,0,400,"Identity",density = bool)
	 ggsave(filename=paste("plots/question1/eig_ident_small_",density,dotplot,".pdf"), 
	 			 plot = eig_ident.ls[[3]], 
	 			 device = cairo_pdf, 
	 			 width = 210, 
	 			 height = 297, 
	 			 units = "mm")
	 ggsave(filename=paste("plots/question1/eig_ident_large_",density,dotplot,".pdf"), 
	 			 plot = eig_ident.ls[[4]], 
	 			 device = cairo_pdf, 
	 			 width = 210, 
	 			 height = 297, 
	 			 units = "mm")
	 # As expected, largest exponential eigenvalue = p, smallest tends towards 0 as p increases
	 eigen_range_cov_structures$exponential$exp_decay
	 eig_linear.ls <- plot_eigens(eigen_range_ml,eigen_range_diag_ml,eigen_range_LedoitWolf,400,800,"Linear",density = bool)
	 ggsave(filename=paste("plots/question1/eig_linear_small_",density,dotplot,".pdf"), 
	 			 plot = eig_linear.ls[[3]], 
	 			 device = cairo_pdf, 
	 			 width = 250, 
	 			 height = 250, 
	 			 units = "mm")
	 ggsave(filename=paste("plots/question1/eig_linear_large_",density,dotplot,".pdf"), 
	 			 plot = eig_ml_linear.ls[[4]], 
	 			 device = cairo_pdf, 
	 			 width = 250, 
	 			 height = 250, 
	 			 units = "mm")
	 # As expected, largest exponential eigenvalue = p, smallest tends towards 0 as p increases
	 eigen_range_cov_structures$exponential$exp_decay
	 eig_expon.ls <- plot_eigens(eigen_range_ml,eigen_range_diag_ml,eigen_range_LedoitWolf,800,1200,"Exponential",density = bool)
	 ggsave(filename=paste("plots/question1/eig_expon_small_",density,dotplot,".pdf"), 
	 			 plot = eig_expon.ls[[3]], 
	 			 device = cairo_pdf, 
	 			 width = 250, 
	 			 height = 250, 
	 			 units = "mm")
	 ggsave(filename=paste("plots/question1/eig_expon_large_",density,dotplot,".pdf"), 
	 			 plot = eig_expon.ls[[4]], 
	 			 device = cairo_pdf, 
	 			 width = 250, 
	 			 height = 250, 
	 			 units = "mm")
	 # we see generally, smallest eigenvalue is ~ 0.33, largest is ~2.9
	 eigen_range_cov_structures$toeplitz$toeplitz
	 eig_toeplitz.ls <- plot_eigens(eigen_range_ml,eigen_range_diag_ml,eigen_range_LedoitWolf,1200,1600,"Toeplitz",density = bool)
	 ggsave(filename=paste("plots/question1/eig_toeplitz_small_",density,dotplot,".pdf"), 
	 			 plot = eig_toeplitz.ls[[3]], 
	 			 device = cairo_pdf, 
	 			 width = 250, 
	 			 height = 250, 
	 			 units = "mm")
	 ggsave(filename=paste("plots/question1/eig_toeplitz_large_",density,dotplot,".pdf"), 
	 			 plot = eig_toeplitz.ls[[4]], 
	 			 device = cairo_pdf, 
	 			 width = 250, 
	 			 height = 250, 
	 			 units = "mm")
}
 
 
 

 
