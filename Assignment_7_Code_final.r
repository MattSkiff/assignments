mu_values = c(-2.0,-0.5,1.0,2.5) #mean values throughout
prior = c(0.3,0.4,0.2,0.1) #Initial Priors

#OBSERVATION 1 = -0.9
proportional_likelihoods = dnorm(-0.9,mu_values,1) #Proportional Likelihood 2 
proportional_likelihoods

marginal = sum(prior*proportional_likelihoods) #Marginal 1
marginal 

posterior_probabilities = prior*proportional_likelihoods/marginal #Posteriors 1
posterior_probabilities
sum(posterior_probabilities = prior*proportional_likelihoods/marginal) #Posterior check 1

#OBSERVATION 2 = 1.4
prior_2 = posterior_probabilities #Priors 2
prior_2

proportional_likelihoods_2 = dnorm(1.4,mu_values,1) #Proportional Likelihood 2
proportional_likelihoods_2

marginal_2 = sum(proportional_likelihoods_2*prior_2) #Marginal 2
marginal_2

posterior_probabilities_2 = proportional_likelihoods_2*prior_2/marginal_2 #Posteriors 2
posterior_probabilities_2
sum(proportional_likelihoods_2*prior_2/marginal_2) #Posterior check 2

#OBSERVATION 3 = -0.2
priors_3 = posterior_probabilities_2 #Priors 3
priors_3

proportional_likelihoods_3 = dnorm(-0.2,mu_values,1) #Proportional Likelihood 3
proportional_likelihoods_3

marginal_3 = sum(priors_3*proportional_likelihoods_3) #Marginal 3
marginal_3

posterior_probabilities_3 = priors_3*proportional_likelihoods_3/marginal_3 #Posteriors 3
posterior_probabilities_3
sum(proportional_likelihoods_3*priors_3/marginal_3) #Posterior check 3

#OBSERVATION 4 = 0.9
priors_4 = posterior_probabilities_3 #Priors 4
priors_4

proportional_likelihoods_4 = dnorm(0.9,mu_values,1) #Proportional Likelihood 4
proportional_likelihoods_4

marginal_4 = sum(proportional_likelihoods_4*priors_4) #Marginal 4
marginal_4

posterior_probabilities_4 = proportional_likelihoods_4*priors_4/marginal_4 #Posteriors 4
posterior_probabilities_4
sum(proportional_likelihoods_4*priors_4/marginal_4) #Posterior check 4

#USING AVERAGE of OBSERVATIONS y1 -> y4

mean_sample = (0.9-0.9-0.2+1.4)/4
variance_sample = 1^2/4 #Variance of sample is prior std dev squared over sample size
std_dev_sample = sqrt(variance_sample) #Standard Deviation of Sample

proportional_likelihoods_means = dnorm(mean_sample,mu_values,std_dev_sample) #Proportional Likelihood Means
proportional_likelihoods_means

#Likelihood data ~ Likelihood of sample mean (which is Normally distributed with u and var = o^2/n)

marginal_means = sum(prior*proportional_likelihoods_means) #Marginal means
marginal_means

posterior_probabilities_means = prior*proportional_likelihoods_means/marginal_means #Posteriors Means
posterior_probabilities_means
sum(proportional_likelihoods_means*prior/marginal_means) #Posterior check Means

#Q1
proportional_likelihoods
prior*proportional_likelihoods #P x PL 1
marginal
posterior_probabilities

#Q2
prior_2
proportional_likelihoods_2
proportional_likelihoods_2*prior_2 #P x PL 2
marginal_2
posterior_probabilities_2

#Q3
priors_3
proportional_likelihoods_3
priors_3*proportional_likelihoods_3 #P x PL 3
marginal_3
posterior_probabilities_3

#Q4
priors_4
proportional_likelihoods_4
proportional_likelihoods_4*priors_4 #P x PL 4
marginal_4
posterior_probabilities_4

#Q5
prior
proportional_likelihoods_means
prior*proportional_likelihoods_means #P x PL Means
marginal_means
posterior_probabilities_means

#Q6
mean_posterior = sum(posterior_probabilities_means*mu_values) #Expectation of the variance of a random variable
variance_posterior = sum(posterior_probabilities*mu_values^2)-mean_posterior^2

mean_posterior
variance_posterior

#Mean and variance of continuous distribution will be the same as mean and variance of discrete distribution
#Therefore Mean(continuous prior normal distribution);
Mean_CPND = sum(mu_values*prior)
Mean_CPND

#Therefore Variance(continuous prior normal distribution);
Var_CPND = sum((mu_values^2)*prior) - Mean_CPND^2
Var_CPND


#Q7
#Using updating rules

mean_posterior_continuous_prior = ((Var_CPND*mean_sample)+(variance_sample*Mean_CPND))/(Var_CPND+variance_sample)
mean_posterior_continuous_prior

variance_posterior_continuous_prior = (Var_CPND*variance_sample)/(variance_sample+Var_CPND)
variance_posterior_continuous_prior








