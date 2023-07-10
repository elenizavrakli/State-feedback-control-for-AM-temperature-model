library(Matrix) 
library(e1071)
library(pracma)
library(devtools)
library(expm)
library(MASS)
library(tidyverse)
library(netcontrol)

#Setting the seed so that the results are reproducible 
set.seed(10)

n = 6 #number of states 
m = 7 #number of inputs

layer_number = 100 #total number of layers simulated / total time steps

gamma = 0.99 #discount factor for perfomance index

# Define system matrices 
# State matrix
A <- matrix(c(0.992,0.0018,0,0,0,0,
              0.0023,0.9919,0.0043,0,0,0,
              0,-0.0042,1.0009,0.0024,0,0,
              0,0,0.0013,0.9979,0,0,
              0,0,0,0,0.9972,0,
              0,0,0,0,0,0.9953),
            nrow = n, byrow = TRUE)
# Input matrix
B <- matrix(c(1.0033,0,0,0,0,0,-0.2175,
              0,1.0460,0,0,0,0,-0.0788,
              0,0,1.0326,0,0,0,-0.0020,
              0,0,0,0.4798,0,0,-0.0669,
              0,0,0,0,0.8882,0,0.1273,
              0,0,0,0,0,1.1699,-0.1792),
            nrow = n, ncol = m, byrow = TRUE)

#State matrix of augmented system
Tmatr<- as.matrix(bdiag(A,diag(n)))

#Input matrix of augmented system
B1<-as.matrix(rbind(B,matrix(rep(0,n*m),n,m)))

#Weighting matrices for performance index
Q <- diag(n)
C1<-as.matrix(cbind(diag(n),-diag(n)))
Q1 <- t(C1)%*%Q%*%C1 
R <- diag(m)

#Initialize the state 
x <- rep(50,n)
xvec <- x #vector storing past states, used to plot the trajectories

#Define optimal state to be tracked
xopt <- c(155,160,165,170,180,190)

#Define augmented state
X<-c(x,xopt) 

# Initialise control policy randomly
i=0
Knew <- matrix(rnorm(m*2*n,0,10),nrow=m)
K <- matrix(rep(3,m*2*n),nrow=m)

#Obtain optimal control policy through repeated solutions of the Lyapunov equation
while(norm(K-Knew)>0.1 ){
  K = Knew
  P <- dlyap(gamma*(Tmatr-B1%*%K),Q1+t(K)%*%R%*%K)
  Knew <- gamma*ginv(R+gamma*t(B1)%*%P%*%B1)%*%t(B1)%*%P%*%Tmatr 
  i=i+1
}

costvec<-NULL #vector to save the cost terms
for(layr in 1:layer_number){
  #Calculate the optimal control
  u <- - Knew%*%X
  
  #Simulation from the state space system, applying the optimal control 
  xnew<-Tmatr%*%X+B1%*%u
  cost <- t(X)%*%Q1%*%X+t(u)%*%R%*%u
  X<-xnew
  
  #Save the new state and cost terms to the respective memories
  xvec<-cbind(xvec,xnew[1:n])
  costvec<-cbind(costvec,cost)
}

#Prepare data for plotting
xvec<-t(xvec)
xvec<-as.data.frame(cbind(1:(layer_number+1),xvec))
colnames(xvec) <- c("time","heater 1","heater 2","heater 3","heater 4","heater 5","heater 6")

#Plot the trajectories using ggplot
xvec[-1,c("time","heater 1","heater 2","heater 3","heater 4","heater 5","heater 6")] %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=1)+
  theme_classic()+
  labs(x='time',y='heater temperature evolution')
