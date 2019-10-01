## Illustration of discrete Markov chains  ##
# transition prob matrix, long run distrn  ##
# - Chaitanya Joshi                        ##
#############################################

#This is a discrete Markov chain with 3 states - 1, 2 and 3.

P1 = c((3/9),0.4,(4/15))    #p11,p12,p13
P2 = c(0.24,0.56,0.2)  #p21,p22,p23
P3 = c(0.4,0.5,0.1)    #p31,p32,p33

P = rbind(P1,P2,P3)    # 1-step transition prob matrix

P2 = P%*%P             # 2-step transition prob matrix
P2

P3 = P2 %*%P          # 3-step transition prob matrix
P3


P5 = P3%*%P2           # 5-step transition prob matrix
P5

n = 10
P_int = P              #P_intermediate for calculations

n_grand = 1000
n_grand.vec <- NULL
g <- 0

for (i in 1:n_grand) {

g = g+1  
n <- n_grand  
  
for (i in 1:n+1)
{	P_int = P%*%P_int
}
Pn = P_int            # n=step transition prob matrix

n_grand.vec[n_grand] <- Pn

}

a0 = c(0.25,0.35,0.4)   #alpha_0 occupation prob at time 0
an = a0%*%Pn            #alpha_n occuptation prob at time n
an

pi = c(0.2992701, 0.3065693, 0.3941606)   #for n>10 'an' remains constant - long run distribution

pi%*%P                #pi = pi*P  is satistfied

mii = 1/pi             #mean recurrrence time for each state
mii

### 

P1 = c(0.75,0.1,0.15)    #p11,p12,p13
P2 = c(0.06,0.9,0.04)  #p21,p22,p23
P3 = c(0.225,0.1,0.675)    #p31,p32,p33

P = rbind(P1,P2,P3)    # 1-step transition prob matrix

P2 = P%*%P             # 2-step transition prob matrix
P2

P3 = P2 %*%P          # 3-step transition prob matrix
P3


P5 = P3%*%P2           # 5-step transition prob matrix
P5

n = 10000
P_int = P              #P_intermediate for calculations

for (i in 1:n)
{	P_int = P%*%P_int
}
Pn = P_int            # n=step transition prob matrix
Pn

