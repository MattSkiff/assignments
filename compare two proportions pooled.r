n1 <- 200 #number of insects in sample A
n2 <- 200 #number of insects in sample B
y1 <- 120 #dead insects in sample A
y2 <- 90 #dead insects in sample B

p1 <- y1 / n1
p2 <- y2 / n2

common_pop_prop <- ((n1*p1)+(n2*p2))/(n1+n2)

z <- (p1-p2) / sqrt(common_pop_prop*(1-common_pop_prop)*((1/n1)+(1/n2)))
z
