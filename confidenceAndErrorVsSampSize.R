# Evaluating error and confidence as functions of the sample size.
# For a 1-\alpha C.I. we usually choose sample size:
# n > z_{\alpha/2}*S^2/\epsilon^2 (*)
# where \epsilon and \alpha are user-defined error and confidence 
# parameters respectively.
# Interesting to see how these two quantities depend on the sample size.
# For convenience in (*) we will set S^2 = 1. For studying error behaviour,
# we set \alpha = 0.025 (i.e. 95% confidence) and for studying confidence
# we set \epsilon = 0.1 and we see how the other term varies as a function
# of sample size.

rm(list=ls()) #clear the environment
plot.new() #create new plot object
par(mfrow = c(2,2)) #split the plot into 4 subplots in a 2*2 array, filled by row (mfcol fills by column)
sampSizes = c(1:1000)

# Error

zCrit = qnorm(0.975)
error = zCrit/sqrt(sampSizes)
plot(sampSizes,error,col='green') #col='green' changes colour from the default black to green
plot(sampSizes,2/sqrt(sampSizes),type='l',col='red',lwd=2) # type='l' makes a line plot, lwd = 2 makes line twice usual thickness


# Confidence
epsilon = 0.1
confidence = pnorm(epsilon*sqrt(sampSizes)) #pnorm is the inverse of qnorm
plot(sampSizes,confidence,col='blue')
plot(sampSizes,1-exp(-sampSizes*epsilon^2),type='l',col='black',lwd=2)

