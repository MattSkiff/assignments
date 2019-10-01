xBarN = numeric(100)
for(i in 1:100)
{
xn = rnorm(100)
# hist(xn)
xBarN[i] = mean(xn)
}
hist(xBarN)

xBarU = numeric(100)
for(i in 1:100)
{
xu = runif(100)
# hist(xu)
xBarU[i] = mean(xu)
}
hist(xBarU)

xBarC2 = numeric(100)
for(i in 1:100)
{
xc2 = rchisq(100, 1, ncp = 0)
# hist(xc2)
xBarC2[i] = mean(xc2)
}
hist(xBarC2)

