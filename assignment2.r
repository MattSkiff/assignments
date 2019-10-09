library(MASS) # Generation of samples from multivariate normal distributions
library(corpcor) # Ledoit Wolf shrinkage estimator
library(nlshrink) # 2nd ver Ledoit Wolf shrinkage estimator
library(ddpcr) # Gets cov.shrink() to STFU
library(purrr) # calculating dot products by mapping sample matrices to covariance
library(ggplot2);library(ggridges) # customised graphics
library(grid);library(gridExtra) # arranging arrays of ggplots
library(reshape2) # melting dataframes
library(stringr) # str_extract
library(lemon) # grid_arrange_shared_legend()
library(naturalsort) #

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
## Here, we use three types of estimators: maximum likelihood with known mean, a diagonal estimator of only variances, and the Ledoit-Wolf Shrinkage shrinkage estimator

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
saveRDS(diag_ml_cov_estimates,"ml_cov_estimates.rds")
saveRDS(ml_cov_estimates,"diag_ml_cov_estimates.rds")
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
 
 # Ledoit-Wolf Shrinkage
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
 
 # Ledoit-Wolf Shrinkage
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
 
# A casual survey of the summary statistics of the estimates against the ground truth does not reveal much; hence a more in depth visuali inspection is needed.
 
 summary(rapply(trace_ml,identity))
 summary(rapply(trace_diag_ml,identity))
 summary(rapply(trace_LedoitWolf,identity))
 summary(rapply(trace_cov_structures,identity))
 
#---- Plotting ----
 
#---- Q1 Plots ----

# defining common ggproto objects to reduce verbosity 
rug_ <- geom_rug(mapping = aes(x = value, colour = variable),show.legend = F,alpha = 1/2) 
ridges_ <- stat_density_ridges(mapping = aes(x = value,y = variable,fill = variable),color = "black",quantile_lines = T,quantiles = 2)
theme_ridge <- theme_light() + theme(axis.text.x = element_text(angle = 0, hjust = 1,  vjust = -0.5,size = 7),axis.text.y = element_blank(),axis.title = element_text(size = 7)) +	theme(plot.title = element_text(size = 7))
# + guides(color = FALSE) + scale_y_discrete(expand = expand_scale(add = c(0.2, 2.5)))
 
 # ML,D-ML,LW Plotting Function

plot_eigens <- function(x,y,z,start,end,covariance_structure,density = F) {
	
	s <- switch (covariance_structure,
							 "Identity" = 1,
							 "Linear" = 2,
							 "Exponential" = 3,
							 "Toeplitz" = 4
	)
	
	ref <- my_flatten(eigen_range_cov_structures[[s]])
	
	options(digits = 3, scipen = -2)
	 eigen_smallest <- list()
	 eigen_biggest <- list()
	 df_small <- data.frame(matrix(NA, ncol=1, nrow=21)[-1]) # from SO
	 df_large <- data.frame(matrix(NA, ncol=1, nrow=21)[-1]) # from SO
	 for (i in 1:20) {
	 	title <- names(my_flatten(x)[start+i])
	 	title <- str_remove(pattern = str_extract(pattern = "[^.]*.[^.]*",string = title),string = title)
	 	title <- substr(title,2,nchar(title))
	 	df_x <- as.data.frame(do.call(rbind, my_flatten(x)[seq(start+i,end,by = 20)]))
	 	df_y <- as.data.frame(do.call(rbind, my_flatten(y)[seq(start+i,end,by = 20)]))
	 	df_z <- as.data.frame(do.call(rbind, my_flatten(z)[seq(start+i,end,by = 20)]))
	 	colnames(df_x) <- c('Smallest_ML',"Largest_ML")
	 	colnames(df_y) <- c('Smallest_D_ML',"Largest_D_ML")
	 	colnames(df_z) <- c('Smallest_TP',"Largest_TP")
	 	quiet(df <- melt(cbind(df_x,df_y,df_z)))
	 	df_small <- df[(df[,1] == "Smallest_ML" | df[,1] == "Smallest_D_ML" | df[,1] == "Smallest_TP" ),]
	 	df_large <- df[(df[,1] == "Largest_ML" | df[,1] == "Largest_D_ML" | df[,1] == "Largest_TP" ),]
	 	v <- as.numeric(ref[[rep(seq(1:5),4)[i]]])
			if (density) {
			
				g <- ggplot(data = df_small) + 
					ridges_+
					rug_ +
					geom_vline(xintercept = v[1],color = 'red',show.legend = F) +
					labs(title = title,y = "density", x = "smallest eigenvalues") +
					theme_ridge + guides(color = FALSE) + scale_y_discrete(expand = expand_scale(add = c(0.2, 2.5))) + # scal y fixes tops of ridges being cutoff
					scale_fill_discrete(name = paste(covariance_structure," Structure |  Estimator: "),labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage"))
				
				gl <- ggplot(data = df_large) + 
					ridges_+
					rug_ +
					geom_vline(xintercept = v[2],color = 'red',show.legend = F) +
					labs(title = title,y = "density", x = "largest eigenvalues") +
					theme_ridge + guides(color = FALSE) + scale_y_discrete(expand = expand_scale(add = c(0.2, 2.5))) +
					scale_fill_discrete(name = paste(covariance_structure," Structure |  Estimator:"),labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage"))
			} else if (!density) {
			if (covariance_structure == "Identity") {
		 		g <- ggplot(data = df_small) + 
		 			geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
		 			geom_vline(xintercept = v[1],color = 'red') +
		 			labs(title = title,y = "count", x = "smallest eigenvalues") +
		 			theme_light() +
		 			theme(axis.text.x = element_text(angle = 0, hjust = 1,  vjust = -0.5,size = 7),axis.text.y = element_text(size = 7),axis.title = element_text(size = 7)) +
		 			guides(color = FALSE) +
		 			theme(plot.title = element_text(size = 7)) + scale_fill_discrete(name = "Identity Structure | Estimator",labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage"))
		 	} else if (covariance_structure == "Exponential") {
			 	g <- ggplot(data = df_small) + 
			 		geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
			 		geom_vline(xintercept = v[1],color = 'red') +
			 		labs(title = title,y = "count", x = "smallest eigenvalues") +
			 		theme_light() +
			 		theme(plot.title = element_text(size = 7)) +
			 		guides(color = FALSE) +
			 		theme(axis.text.x = element_text(angle = 0, hjust = 1,  vjust = -0.5,size = 7),axis.text.y = element_text(size = 7),axis.title = element_text(size = 7)) + 
					scale_fill_discrete(name = "Exponential Structure | Estimator",labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage"))
		 	} else {
		 		g <- ggplot(data = df_small) + 
		 			geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
		 			geom_vline(xintercept = v[1],color = 'red') +
		 			labs(title = title,y = "count", x = "smallest eigenvalues") +
		 			theme_light() +
		 			theme(axis.text.x = element_text(angle = 0, hjust = 1,  vjust = -0.5,size = 7),axis.text.y = element_text(size = 7),axis.title = element_text(size = 7)) +
		 			guides(color = FALSE) +
		 			theme(plot.title = element_text(size = 7)) + scale_fill_discrete(name = paste(covariance_structure," Structure |  Estimator:"),labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage"))
		 	}
		 	
		 	if (covariance_structure == "Identity") {
		 		gl <- ggplot(data = df_large) + 
		 			geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
		 			geom_vline(xintercept = v[2],color = 'red') +
		 			labs(title = title,y = "count", x = "largest eigenvalues") +
					theme_light() +
		 			guides(color = FALSE) +
		 			theme(axis.text.x = element_text(angle = 0, hjust = 1,  vjust = -0.5,size = 7),axis.text.y = element_text(size = 7),axis.title = element_text(size = 7)) +
		 			theme(plot.title = element_text(size = 7)) + 
					scale_fill_discrete(name = "Identity Structure | Estimator",labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage"))	# a mad daft hack: put the title in the legend
		 	} else if (covariance_structure == "Exponential") {
		 		gl <- ggplot(data = df_large) + 
		 			geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
		 			geom_vline(xintercept = v[2],color = 'red') +
		 			labs(title = title,y = "count", x = "largest eigenvalues") +
		 			theme_light() +
		 			guides(color = FALSE) +
		 			theme(plot.title = element_text(size = 7)) +
		 			theme(axis.text.x = element_text(angle = 0, hjust = 1,  vjust = -0.5,size = 7),axis.text.y = element_text(size = 7),axis.title = element_text(size = 7)) +
					scale_fill_discrete(name = "Exponential Structure | Estimator: ",labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage"))
					
		 	} else {
		 		gl <- ggplot(data = df_large) + 
		 			geom_dotplot(mapping = aes(x = value,color = variable,fill = variable),stackgroups = T,method = "histodot") +
		 			geom_vline(xintercept = v[2],color = 'red') +
		 			labs(title = title,y = "count", x = "largest eigenvalues") +
		 			theme_light() +
		 			guides(color = FALSE) +
		 			theme(axis.text.x = element_text(angle = 0, hjust = 1,  vjust = -0.5,size = 7),axis.text.y = element_text(size = 7),axis.title = element_text(size = 7)) +
		 			theme(plot.title = element_text(size = 7)) + scale_fill_discrete(name = paste(covariance_structure," Structure |  Estimator: "),labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage")) 
			 	}
			} 
	 	eigen_smallest[[i]] <- g
	 	eigen_biggest[[i]] <- gl
	 }
	 density_ <- NULL;dotplot_ <- NULL
	 if (density) { 
	 	density_ <- "KDEs"
	 	dotplot_ <- ""
	 } else { 
	 	dotplot_ <- "Dotplots" 
	 	density_ <- ""
	 }
	 
	 nplots = 20
	 ncol = 5
	 nrow = 4
	 eval(parse(text = paste0("quiet(grob <- grid_arrange_shared_legend(", paste0("eigen_smallest", "[[", c(1:nplots), "]]", sep = '', collapse = ','), ",ncol =", ncol, ",nrow =", nrow, ", position = 'bottom',  
	 												 top=grid::textGrob(paste(dotplot_,density_,'of Smallest Eigenvalues from ',covariance_structure,' Covariance Estimates'), gp=grid::gpar(fontsize=12)),bottom=grid::textGrob('Red lines indicate true population covariance smallest eigenvalues',gp = gpar(fontface = 3, fontsize = 9),hjust = 1,x = 1),plot = F))", sep = '')))
	 eval(parse(text = paste0("quiet(grob_large <- grid_arrange_shared_legend(", paste0("eigen_biggest", "[[", c(1:nplots), "]]", sep = '', collapse = ','), ",ncol =", ncol, ",nrow =", nrow, ", position = 'bottom',  
	 												 top=grid::textGrob(paste(dotplot_,density_,'of Largest Eigenvalues from ',covariance_structure,' Covariance Estimates'), gp=grid::gpar(fontsize=12)),bottom=grid::textGrob('Red lines indicate true population covariance largest eigenvalues',gp = gpar(fontface = 3, fontsize = 9),hjust = 1,x = 1),plot = F))", sep = '')))
	 

	 return(list('small' = df_small,'large' = df_large,'g' = grob,'gl' = grob_large))
}
	# we expect largest and smallest values to be 1
direct <- NULL
for (i in 1:2) {
	if (i == 1) { 
		bool = F
		density <- ""
		dotplot <- "dotplots"
		direct <- "density_plots/"
	}	else {
		 bool = T
		 density <- "density"
		 dotplot <- ""
		 direct <- "density_plots/"
	}
	 eig_ident.ls <- plot_eigens(eigen_range_ml,eigen_range_diag_ml,eigen_range_LedoitWolf,0,400,"Identity",density = bool)
	 ggsave(filename=paste("plots/question1/",direct,"eig_ident_small_",density,dotplot,".pdf"), 
	 			 plot = eig_ident.ls[[3]], 
	 			 device = cairo_pdf, 
	 			 width = 250, 
	 			 height = 250, 
	 			 units = "mm")
	 ggsave(filename=paste("plots/question1/",direct,"eig_ident_large_",density,dotplot,".pdf"), 
	 			 plot = eig_ident.ls[[4]], 
	 			 device = cairo_pdf, 
	 			 width = 250, 
	 			 height = 250, 
	 			 units = "mm")
	 # As expected, largest exponential eigenvalue = p, smallest tends towards 0 as p increases
	 eigen_range_cov_structures$exponential$exp_decay
	 eig_linear.ls <- plot_eigens(eigen_range_ml,eigen_range_diag_ml,eigen_range_LedoitWolf,400,800,"Linear",density = bool)
	 ggsave(filename=paste("plots/question1/",direct,"eig_linear_small_",density,dotplot,".pdf"),
	 			 plot = eig_linear.ls[[3]],
	 			 device = cairo_pdf,
	 			 width = 250,
	 			 height = 250,
	 			 units = "mm")
	 ggsave(filename=paste("plots/question1/",direct,"eig_linear_large_",density,dotplot,".pdf"),
	 			 plot = eig_linear.ls[[4]],
	 			 device = cairo_pdf,
	 			 width = 250,
	 			 height = 250,
	 			 units = "mm")
	 # As expected, largest exponential eigenvalue = p, smallest tends towards 0 as p increases
	 eigen_range_cov_structures$exponential$exp_decay
	 eig_expon.ls <- plot_eigens(eigen_range_ml,eigen_range_diag_ml,eigen_range_LedoitWolf,800,1200,"Exponential",density = bool)
	 ggsave(filename=paste("plots/question1/",direct,"eig_expon_small_",density,dotplot,".pdf"),
	 			 plot = eig_expon.ls[[3]],
	 			 device = cairo_pdf,
	 			 width = 250,
	 			 height = 250,
	 			 units = "mm")
	 ggsave(filename=paste("plots/question1/",direct,"eig_expon_large_",density,dotplot,".pdf"),
	 			 plot = eig_expon.ls[[4]],
	 			 device = cairo_pdf,
	 			 width = 250,
	 			 height = 250,
	 			 units = "mm")
	 # we see generally, smallest eigenvalue is ~ 0.33, largest is ~2.9
	 eigen_range_cov_structures$toeplitz$toeplitz
	 eig_toeplitz.ls <- plot_eigens(eigen_range_ml,eigen_range_diag_ml,eigen_range_LedoitWolf,1200,1600,"Toeplitz",density = bool)
	 ggsave(filename=paste("plots/question1/",direct,"eig_toeplitz_small_",density,dotplot,".pdf"),
	 			 plot = eig_toeplitz.ls[[3]],
	 			 device = cairo_pdf,
	 			 width = 250,
	 			 height = 250,
	 			 units = "mm")
	 ggsave(filename=paste("plots/question1/",direct,"eig_toeplitz_large_",density,dotplot,".pdf"),
	 			 plot = eig_toeplitz.ls[[4]],
	 			 device = cairo_pdf,
	 			 width = 250,
	 			 height = 250,
	 			 units = "mm")
}

#---- Q2 Plots ----

# As some of the kernel density estimates overlapped nearly perfectly, decided to use ridgeline plots to clearly show all three estimators

plot_traces <- function(x,y,z,start,end,covariance_structure) {
	options(digits = 3, scipen = -2)

	df_trace <- data.frame(matrix(NA, ncol=1, nrow=21)[-1]) # from SO
	
	s <- switch (covariance_structure,
					"Identity" = 1,
					"Linear" = 2,
					"Exponential" = 3,
					"Toeplitz" = 4
	)
	
	ref <- my_flatten(trace_cov_structures[[s]])
	traces <- list()
	
	for (i in 1:20) {
		title <- names(my_flatten(x)[start+i])
		title <- str_remove(pattern = str_extract(pattern = "[^.]*.[^.]*",string = title),string = title)
		title <- substr(title,2,nchar(title))
		df_x <- as.data.frame(do.call(rbind, my_flatten(x)[seq(start+i,end,by = 20)]))
		df_y <- as.data.frame(do.call(rbind, my_flatten(y)[seq(start+i,end,by = 20)]))
		df_z <- as.data.frame(do.call(rbind, my_flatten(z)[seq(start+i,end,by = 20)]))
		

		colnames(df_x) <- 'Trace_ML'
		colnames(df_y) <- 'Trace_D_ML'
		colnames(df_z) <- 'Trace_TP'
		v <- as.numeric(ref[[rep(seq(1:5),4)[i]]])
		quiet(df_trace <- melt(cbind(df_x,df_y,df_z)))
			g <- ggplot(data = df_trace) + 
				ridges_+
				geom_vline(xintercept = v,colour = "red",show.legend = F) +
				rug_ +
				labs(title = title,y = "density", x = "trace") +
				theme_ridge + + guides(color = FALSE) + scale_y_discrete(expand = expand_scale(add = c(0.2, 2.5))) +
				scale_fill_discrete(name = "Estimator:",labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage"))
		traces[[i]] <- g
	}

# It's digusting, but it works...
# https://stackoverflow.com/questions/6496811/how-to-pass-a-list-to-a-function-in-r (comments)
nplots = 20
ncol = 5
nrow = 4
eval(parse(text = paste0("quiet(grob <- grid_arrange_shared_legend(", paste0("traces", "[[", c(1:nplots), "]]", sep = '', collapse = ','), ",ncol =", ncol, ",nrow =", nrow, ", position = 'bottom',  
													top=grid::textGrob(paste('KDEs of Population Trace',covariance_structure,' Covariance Estimates'), gp=grid::gpar(fontsize=12)),bottom=grid::textGrob('Red lines indicate true population covariance trace',gp = gpar(fontface = 3, fontsize = 9),hjust = 1,x = 1),plot = F))", sep = '')))
	
return(list('trace' = df_trace,'g' = grob))
}

# Trace of identity matrix = p
trace_cov_structures$identity$identity
trace_ident.ls <- plot_traces(trace_ml,trace_diag_ml,trace_LedoitWolf,0,400,"Identity")
ggsave(filename=paste("plots/question2/trace_ident_density.pdf",density,dotplot,".pdf"), 
			 plot = trace_ident.ls[[2]], 
			 device = cairo_pdf, 
			 width = 250, 
			 height = 250, 
			 units = "mm")
# Trace of linear will simply be sum(1:p)
trace_cov_structures$linear$linear_decay
trace_linear.ls <- plot_traces(trace_ml,trace_diag_ml,trace_LedoitWolf,400,800,"Linear")
ggsave(filename=paste("plots/question2/trace_linear_density.pdf",density,dotplot,".pdf"),
			 plot = trace_linear.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")
# As expected, largest exponential eigenvalue = p, smallest tends towards 0 as p increases
trace_cov_structures$exponential$exp_decay
trace_expon.ls <- plot_traces(trace_ml,trace_diag_ml,trace_LedoitWolf,800,1200,"Exponential")
ggsave(filename=paste("plots/question2/trace_expon_density.pdf"),
			 plot = trace_expon.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")
# we see generally, smallest eigenvalue is ~ 0.33, largest is ~2.9
trace_cov_structures$toeplitz$toeplitz
trace_toeplitz.ls <- plot_traces(trace_ml,trace_diag_ml,trace_LedoitWolf,1200,1600,"Toeplitz")
ggsave(filename=paste("plots/question2/trace_toeplitz_density.pdf",density,dotplot,".pdf"),
			 plot = trace_toeplitz.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")

#---- Q3 Plots ----

# boxplots
options(digits = 5, scipen = -2)
# As some of the kernel density estimates overlapped nearly perfectly, decided to use ridgeline plots to clearly show all three estimators
dp_diag_ml <- list('identity' = dot_products_diag_ml_identity,'linear' = dot_products_diag_ml_linear,'exponential' = dot_products_diag_ml_exponential,'toeplitz' = dot_products_diag_ml_toeplitz)
dp_ml <- list('identity' = dot_products_ml_identity,'linear' = dot_products_ml_linear,'exponential' = dot_products_ml_exponential,'toeplitz' = dot_products_ml_toeplitz)
dp_LedoitWolf <- list('identity' = dot_products_LedoitWolf_cov_estimates_identity,'linear' = dot_products_LedoitWolf_cov_estimates_linear,'exponential' = dot_products_LedoitWolf_cov_estimates_exponential,'toeplitz' = dot_products_LedoitWolf_cov_estimates_toeplitz)

plot_dps <- function(x,y,z,start,end,covariance_structure) {
	options(digits = 3, scipen = -2)
	# x is estimator list of eigens
	dps <- list()
	
	df_dp <- data.frame(matrix(NA, ncol=1, nrow=21)[-1]) # from SO
	
	s <- switch (covariance_structure,
							 "Identity" = 1,
							 "Linear" = 2,
							 "Exponential" = 3,
							 "Toeplitz" = 4
	)
	
	dps <- list()
	
	for (i in 1:20) {
		title <- names(my_flatten(x)[start+i])
		title <- str_remove(pattern = str_extract(pattern = "[^.]*.[^.]*",string = title),string = title)
		title <- substr(title,2,nchar(title))
		df_x <- as.data.frame(do.call(rbind, my_flatten(x)[seq(start+i,end,by = 20)]))
		df_y <- as.data.frame(do.call(rbind, my_flatten(y)[seq(start+i,end,by = 20)]))
		df_z <- as.data.frame(do.call(rbind, my_flatten(z)[seq(start+i,end,by = 20)]))
		
		
		colnames(df_x) <- 'dp_ML'
		colnames(df_y) <- 'dp_D_ML'
		colnames(df_z) <- 'dp_TP'
		quiet(df_dp <- melt(cbind(df_x,df_y,df_z)))
		g <- ggplot(data = df_dp) + 
			geom_boxplot(mapping = aes(y = value,x = variable,fill = variable),color = alpha("black", 0.1),outlier.shape = NA,alpha = 0.5) +
			geom_hline(yintercept = 1,color = 'red') +
			geom_jitter(mapping = aes(y = value,x = variable,fill = variable,color = variable),width = 0.4,size = 0.5, alpha = 0.8) +
			geom_rug(mapping = aes(x = value, colour = variable),show.legend = F,alpha = 1/2,sides = 'r') +
			labs(title = title,y = "dot products", x = "estimator") +
			theme_light() + theme(axis.text.y = element_text(angle = 0, hjust = 1,  vjust = -0.5,size = 7),axis.text.x = element_blank(),axis.title = element_text(size = 7)) +	theme(plot.title = element_text(size = 7)) + guides(color = FALSE) +
			scale_fill_discrete(name = "Estimator:",labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage")) + ylim(-1,1)
		dps[[i]] <- g
	}
	
	# It's digusting, but it works...
	# https://stackoverflow.com/questions/6496811/how-to-pass-a-list-to-a-function-in-r (comments)
	nplots = 20
	ncol = 5
	nrow = 4
	eval(parse(text = paste0("quiet(grob <- grid_arrange_shared_legend(", paste0("dps", "[[", c(1:nplots), "]]", sep = '', collapse = ','), ",ncol =", ncol, ",nrow =", nrow, ", position = 'bottom',  
													top=grid::textGrob(paste('Boxplots of Dot Products between Principal Eigenvector Estimate and True Eigenvector Estimate from ',covariance_structure, ' Covariance Structure'), gp=grid::gpar(fontsize=10)),plot = F))", sep = '')))
	
	return(list('dp' = df_dp,'g' = grob))
}

dp_ident.ls <- plot_dps(dp_diag_ml,dp_diag_ml,dp_LedoitWolf,0,400,"Identity")
ggsave(filename=paste("plots/question3/boxplots/dp_ident_boxplots.pdf"), 
			 plot = dp_ident.ls[[2]], 
			 device = cairo_pdf, 
			 width = 250, 
			 height = 250, 
			 units = "mm")

dp_linear.ls <- plot_dps(dp_ml,dp_diag_ml,dp_LedoitWolf,400,800,"Linear")
ggsave(filename=paste("plots/question3/boxplots/dp_linear_boxplots.pdf"),
			 plot = dp_linear.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")

dp_expon.ls <- plot_dps(dp_ml,dp_diag_ml,dp_LedoitWolf,800,1200,"Exponential")
ggsave(filename=paste("plots/question3/boxplots/dp_expon_boxplots.pdf"),
			 plot = dp_expon.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")

dp_toeplitz.ls <- plot_dps(dp_ml,dp_diag_ml,dp_LedoitWolf,1200,1600,"Toeplitz")
ggsave(filename=paste("plots/question3/boxplots/dp_toeplitz_boxplots.pdf"),
			 plot = dp_toeplitz.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")

# density plots
# As some of the kernel density estimates overlapped nearly perfectly, decided to use ridgeline plots to clearly show all three estimators
dp_diag_ml <- list('identity' = dot_products_diag_ml_identity,'linear' = dot_products_diag_ml_linear,'exponential' = dot_products_diag_ml_exponential,'toeplitz' = dot_products_diag_ml_toeplitz)
dp_ml <- list('identity' = dot_products_ml_identity,'linear' = dot_products_ml_linear,'exponential' = dot_products_ml_exponential,'toeplitz' = dot_products_ml_toeplitz)
dp_LedoitWolf <- list('identity' = dot_products_LedoitWolf_cov_estimates_identity,'linear' = dot_products_LedoitWolf_cov_estimates_linear,'exponential' = dot_products_LedoitWolf_cov_estimates_exponential,'toeplitz' = dot_products_LedoitWolf_cov_estimates_toeplitz)

plot_dps <- function(x,y,z,start,end,covariance_structure) {
	options(digits = 3, scipen = -2)
	# x is estimator list of eigens
	dps <- list()
	
	df_dp <- data.frame(matrix(NA, ncol=1, nrow=21)[-1]) # from SO
	
	s <- switch (covariance_structure,
							 "Identity" = 1,
							 "Linear" = 2,
							 "Exponential" = 3,
							 "Toeplitz" = 4
	)
	
	dps <- list()
	
	for (i in 1:20) {
		title <- names(my_flatten(x)[start+i])
		title <- str_remove(pattern = str_extract(pattern = "[^.]*.[^.]*",string = title),string = title)
		title <- substr(title,2,nchar(title))
		df_x <- as.data.frame(do.call(rbind, my_flatten(x)[seq(start+i,end,by = 20)]))
		df_y <- as.data.frame(do.call(rbind, my_flatten(y)[seq(start+i,end,by = 20)]))
		df_z <- as.data.frame(do.call(rbind, my_flatten(z)[seq(start+i,end,by = 20)]))
		
		
		colnames(df_x) <- 'dp_ML'
		colnames(df_y) <- 'dp_D_ML'
		colnames(df_z) <- 'dp_TP'
		quiet(df_dp <- melt(cbind(df_x,df_y,df_z)))
		g <- ggplot(data = df_dp) + 
			stat_density_ridges(mapping = aes(x = value,y = variable,fill = variable),color = "black",quantile_lines = T,quantiles = 2) +
			geom_vline(xintercept = 1,color = 'red') +
			rug_ +
			labs(title = title,y = "density", x = "dot products") +
			theme_ridge + guides(color = FALSE) + scale_y_discrete(expand = expand_scale(add = c(0.2, 2.5))) +
			scale_fill_discrete(name = "Estimator:",labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage"))
		dps[[i]] <- g
	}
	
	# It's digusting, but it works...
	# https://stackoverflow.com/questions/6496811/how-to-pass-a-list-to-a-function-in-r (comments)
	nplots = 20
	ncol = 5
	nrow = 4
	eval(parse(text = paste0("quiet(grob <- grid_arrange_shared_legend(", paste0("dps", "[[", c(1:nplots), "]]", sep = '', collapse = ','), ",ncol =", ncol, ",nrow =", nrow, ", position = 'bottom',  
													top=grid::textGrob(paste('KDEs of Dot Products between Principal Eigenvector Estimate and True Eigenvector Estimate from ',covariance_structure, ' Covariance Structure'), gp=grid::gpar(fontsize=10)),bottom=grid::textGrob('Red lines indicate true population covariance smallest eigenvalues',gp = gpar(fontface = 3, fontsize = 9),hjust = 1,x = 1),plot = F))", sep = '')))
	
	return(list('dp' = df_dp,'g' = grob))
}

dp_ident.ls <- plot_dps(dp_diag_ml,dp_diag_ml,dp_LedoitWolf,0,400,"Identity")
ggsave(filename=paste("plots/question3/density_plots/dp_ident_density.pdf",density,dotplot,".pdf"), 
			 plot = dp_ident.ls[[2]], 
			 device = cairo_pdf, 
			 width = 250, 
			 height = 250, 
			 units = "mm")

dp_linear.ls <- plot_dps(dp_ml,dp_diag_ml,dp_LedoitWolf,400,800,"Linear")
ggsave(filename=paste("plots/question3/density_plots/dp_linear_density.pdf",density,dotplot,".pdf"),
			 plot = dp_linear.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")

dp_expon.ls <- plot_dps(dp_ml,dp_diag_ml,dp_LedoitWolf,800,1200,"Exponential")
ggsave(filename=paste("plots/question3/density_plots/dp_expon_density.pdf"),
			 plot = dp_expon.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")

dp_toeplitz.ls <- plot_dps(dp_ml,dp_diag_ml,dp_LedoitWolf,1200,1600,"Toeplitz")
ggsave(filename=paste("plots/question3/density_plots/dp_toeplitz_density.pdf",density,dotplot,".pdf"),
			 plot = dp_toeplitz.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")

#---- Q4 Plots ----
 
# As some of the kernel density estimates overlapped nearly perfectly, decided to use ridgeline plots to clearly show all three estimators
sse_diag_ml <- list('identity' = SSE_diag_ml_identity,'linear' = SSE_diag_ml_linear,'exponential' = SSE_diag_ml_exponential,'toeplitz' = SSE_diag_ml_toeplitz)
sse_ml <- list('identity' = SSE_ml_identity,'linear' = SSE_ml_linear,'exponential' = SSE_ml_exponential,'toeplitz' = SSE_ml_toeplitz)
sse_LedoitWolf <- list('identity' = SSE_LedoitWolf_cov_estimates_identity,'linear' = SSE_LedoitWolf_cov_estimates_linear,'exponential' = SSE_LedoitWolf_cov_estimates_exponential,'toeplitz' = SSE_LedoitWolf_cov_estimates_toeplitz)

plot_sses <- function(x,y,z,start,end,covariance_structure) {
	options(digits = 3, scipen = -2)
	# x is estimator list of eigens
	sses <- list()
	
	df_sse <- data.frame(matrix(NA, ncol=1, nrow=21)[-1]) # from SO
	
	s <- switch (covariance_structure,
							 "Identity" = 1,
							 "Linear" = 2,
							 "Exponential" = 3,
							 "Toeplitz" = 4
	)
	
	sses <- list()
	
	for (i in 1:20) {
		title <- names(my_flatten(x)[start+i])
		title <- str_remove(pattern = str_extract(pattern = "[^.]*.[^.]*",string = title),string = title)
		title <- substr(title,2,nchar(title))
		#title <- substr(title,17,title) # trim off "identity.identity"
		df_x <- as.data.frame(do.call(rbind, my_flatten(x)[seq(start+i,end,by = 20)]))
		df_y <- as.data.frame(do.call(rbind, my_flatten(y)[seq(start+i,end,by = 20)]))
		df_z <- as.data.frame(do.call(rbind, my_flatten(z)[seq(start+i,end,by = 20)]))
		
		
		colnames(df_x) <- 'sse_ML'
		colnames(df_y) <- 'sse_D_ML'
		colnames(df_z) <- 'sse_TP'
		quiet(df_sse <- melt(cbind(df_x,df_y,df_z)))
		g <- ggplot(data = df_sse) + 
			#geom_density(mapping = aes(x = value,color = variable,fill = variable),alpha = 0.5) +
			#geom_density_ridges2(mapping = aes(x = value,y = variable,color = variable,fill = variable)) +
			ridges_+
			rug_ +
			labs(title = title,y = "density", x = "sse") +
			theme_ridge + guides(color = FALSE) + scale_y_discrete(expand = expand_scale(add = c(0.2, 2.5))) +
			scale_fill_discrete(name = "Estimator:",labels = c("Maximum Likelihood", "Diagonalised Maximum Likelihood", "Ledoit-Wolf Shrinkage"))
		sses[[i]] <- g
	}

	# It's digusting, but it works...
	# https://stackoverflow.com/questions/6496811/how-to-pass-a-list-to-a-function-in-r (comments)
	nplots = 20
	ncol = 5
	nrow = 4
	eval(parse(text = paste0("quiet(grob <- grid_arrange_shared_legend(", paste0("sses", "[[", c(1:nplots), "]]", sep = '', collapse = ','), ",ncol =", ncol, ",nrow =", nrow, ", position = 'bottom',  
													top=grid::textGrob(paste('KDEs of Covariance Estimate SSE from',covariance_structure,' Covariance Structure'), gp=grid::gpar(fontsize=12)),plot = F))", sep = '')))
	
	#quiet(grob <- grid_arrange_shared_legend(sses[[1]],sses[[1]],top = "hello")) # test success
	
	#quiet(grob <- grid_arrange_shared_legend_plotlist(plotlist = sses,ncol = 5,top = "hello"))#,bottom = b_)) #grid.arrange(grobs = sses,ncol = 5) ,top = "test"
	
	return(list('sse' = df_sse,'g' = grob))
}

sse_ident.ls <- plot_sses(sse_ml,sse_diag_ml,sse_LedoitWolf,0,400,"Identity")
ggsave(filename=paste("plots/question4/sse_ident_density.pdf",density,dotplot,".pdf"), 
			 plot = sse_ident.ls[[2]], 
			 device = cairo_pdf, 
			 width = 250, 
			 height = 250, 
			 units = "mm")

sse_linear.ls <- plot_sses(sse_ml,sse_diag_ml,sse_LedoitWolf,400,800,"Linear")
ggsave(filename=paste("plots/question4/sse_linear_density.pdf",density,dotplot,".pdf"),
			 plot = sse_linear.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")

sse_expon.ls <- plot_sses(sse_ml,sse_diag_ml,sse_LedoitWolf,800,1200,"Exponential")
ggsave(filename=paste("plots/question4/sse_expon_density.pdf"),
			 plot = sse_expon.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")

sse_toeplitz.ls <- plot_sses(sse_ml,sse_diag_ml,sse_LedoitWolf,1200,1600,"Toeplitz")
ggsave(filename=paste("plots/question4/sse_toeplitz_density.pdf",density,dotplot,".pdf"),
			 plot = sse_toeplitz.ls[[2]],
			 device = cairo_pdf,
			 width = 250,
			 height = 250,
			 units = "mm")

#---- Summary Plots ---- 


n_sort_list_summary <- function(x) { return(unlist(x)[naturalorder(names(unlist(x)))]) }

unique(names(unlist(dot_products_LedoitWolf_cov_estimates_toeplitz)[naturalorder(names(unlist(dot_products_LedoitWolf_cov_estimates_toeplitz)))]))

# Q3
x11()
pdf("Dot Product Summary.pdf")
par(mfrow = c(3,4))
plot(n_sort_list_summary(dot_products_ml_identity),bg = 'red',main = "Identity - ML",pch = 21,ylab = "Dot Products - ML")
plot(n_sort_list_summary(dot_products_ml_linear),bg = 'blue',main = "Linear - ML",pch = 21,ylab = "Dot Products - ML")
plot(n_sort_list_summary(dot_products_ml_exponential),bg = 'purple',main = "Exponential - ML",pch = 21,ylab = "Dot Products - ML")
plot(n_sort_list_summary(dot_products_ml_toeplitz),bg = 'green',main = "Toeplitz - ML",pch = 21,ylab = "Dot Products - ML")
plot(n_sort_list_summary(dot_products_diag_ml_identity),bg = 'red',main = "Identity - Diag ML",pch = 21,ylab = "Dot Products - Diag ML")
plot(n_sort_list_summary(dot_products_diag_ml_linear),bg = 'blue',main = "Linear - Diag ML",pch = 21,ylab = "Dot Products - Diag ML")
plot(n_sort_list_summary(dot_products_diag_ml_exponential),bg = 'purple',main = "Exponential - Diag ML",pch = 21,ylab = "Dot Products - Diag ML")
plot(n_sort_list_summary(dot_products_diag_ml_toeplitz),bg = 'green',main = "Toeplitz - Diag ML",pch = 21,ylab = "Dot Products - Diag ML")
plot(n_sort_list_summary(dot_products_LedoitWolf_cov_estimates_identity),bg = 'red',main = "Identity - LW",pch = 21,ylab = "Dot Products - LW")
plot(n_sort_list_summary(dot_products_LedoitWolf_cov_estimates_linear),bg = 'blue',main = "Linear - LW",pch = 21,ylab = "Dot Products  - LW")
plot(n_sort_list_summary(dot_products_LedoitWolf_cov_estimates_exponential),bg = 'purple',main = "Exponential - LW",pch = 21,ylab = "Dot Products  - LW")
plot(n_sort_list_summary(dot_products_LedoitWolf_cov_estimates_toeplitz),bg = 'green',main = "Toeplitz - LW",pch = 21,ylab = "Dot Products  - LW")  
dev.off()

# Q4
x11()
pdf("SSE Summary.pdf")
par(mfrow = c(3,4))
plot(n_sort_list_summary(SSE_ml_identity),bg = 'red',main = "Identity - ML",pch = 21,ylab = "SSE - ML")
plot(n_sort_list_summary(SSE_ml_linear),bg = 'blue',main = "Linear - ML",pch = 21,ylab = "SSEs - ML")
plot(n_sort_list_summary(SSE_ml_exponential),bg = 'purple',main = "Exponential - ML",pch = 21,ylab = "SSEs - ML")
plot(n_sort_list_summary(SSE_ml_toeplitz),bg = 'green',main = "Toeplitz - ML",pch = 21,ylab = "SSEs - ML")
plot(n_sort_list_summary(SSE_diag_ml_identity),bg = 'red',main = "Identity - Diag ML",pch = 21,ylab = "SSEs - Diag ML")
plot(n_sort_list_summary(SSE_diag_ml_linear),bg = 'blue',main = "Linear - Diag ML",pch = 21,ylab = "SSEs - Diag ML")
plot(n_sort_list_summary(SSE_diag_ml_exponential),bg = 'purple',main = "Exponential - Diag ML",pch = 21,ylab = "SSEs - Diag ML")
plot(n_sort_list_summary(SSE_diag_ml_toeplitz),bg = 'green',main = "Toeplitz - Diag ML",pch = 21,ylab = "SSEs - Diag ML")
plot(n_sort_list_summary(SSE_LedoitWolf_cov_estimates_identity),bg = 'red',main = "Identity - LW",pch = 21,ylab = "SSEs - LW")
plot(n_sort_list_summary(SSE_LedoitWolf_cov_estimates_linear),bg = 'blue',main = "Linear - LW",pch = 21,ylab = "SSEs  - LW")
plot(n_sort_list_summary(SSE_LedoitWolf_cov_estimates_exponential),bg = 'purple',main = "Exponential - LW",pch = 21,ylab = "SSEs  - LW")
plot(n_sort_list_summary(SSE_LedoitWolf_cov_estimates_toeplitz),bg = 'green',main = "Toeplitz - LW",pch = 21,ylab = "SSEs  - LW")  
dev.off()

#---- Numerical Summaries ----

cbind(c("Identity","Linear","Exponential","Toeplitz"),
			rbind(summary(unlist(dot_products_ml_identity)),
						summary(unlist(dot_products_ml_linear)),
						summary(unlist(dot_products_ml_exponential)),
						summary(unlist(dot_products_ml_toeplitz))))

 
