## Clear the environment:

rm(list=ls())

rows = 10
cols = 10000
mu = 0
sigma = 30
sampPaths = matrix(0,rows,cols)
#head(sampPaths)
sampData = matrix( rnorm(rows*cols,mean=mu,sd=sigma), rows, cols)
#head(sampData)
counter = c(1:cols)
colours <- c(rainbow(rows, s = 1, v = 1, start = 0, end = max(1, rows - 1)/rows, alpha = 1))
for(m in 1:rows)
{
for(n in 1:cols)
{ 
sampPaths[m,n] = mean(sampData[m,1:n])
}
}
plot.new()
plot(counter,sampPaths[1,],type="l",col=colours[1],ann=FALSE, ylim=c(-10, 10))
for(m in 2:rows)
{
lines(counter,sampPaths[m,],col=colours[m])
}
abline(0,0)
title(xlab= "Sample Size")
title(ylab= "Sample Means")

#meanFun <- function(x,start,end){return(mean(x[start:end]))}