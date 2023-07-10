library(Matrix) 
library(e1071)
library(pracma)
library(devtools)
library(expm)
library(tidyverse)
library(MASS)

set.seed(100)

n = 6 #number of states 
m = 7 #number of inputs

layer_number = 100 #total number of layers simulated / total time steps

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

#Weighting matrices for performance index
Q <- diag(n) #weight associated with tracking
R <- diag(m) #weight associated with input

#Initialise memories to store S and b matrices, used for the optimal control calculations
S <- Q #initialise S for the final time point
b = rep(0,n) #initialise b for the final time point

memS <- S #memory storing values of S
memb <- b #memory storing values of b

#Initialize the state 
x <- rep(50,n)
xvec <- x #vector storing past states, used to plot the trajectories

# Define optimal state to be tracked
xopt <- c(155,160,165,170,180,190)

# Solve the equations for b and S backwards in time and store them in memories memb and memS
for(i in 100:1){
  K <- t(-ginv(t(B)%*%S%*%B+R)%*%t(B)%*%S%*%A)
  bnew <- (t(A)+K%*%t(B))%*%(b-Q%*%xopt)
  Snew <- t(A)%*%(S-S%*%B%*%ginv(t(B)%*%S%*%B+R)%*%t(B)%*%S)%*%A+Q
  
  b <- bnew
  S <- Snew
  
  memS <- rbind(S,memS) 
  memb <- c(b,memb)
}

# Simulate the system behaviour using the S(i+1) and b(i+1) in each step, as per the finite horizon algorithm
for(i in 1:98){
  # Retrieve the S and b from the respective memories
  S <- memS[((i)*n+1):((i+1)*n),]
  b <- memb[((i)*n+1):((i+1)*n)]
  
  #Calculate the optimal control
  u <- - ginv(t(B)%*%S%*%B+R)%*%t(B)%*%(S%*%A%*%x+(b-Q%*%xopt))
  
  #Simulation from the state space system, applying the optimal control 
  xnew<-A%*%x+B%*%u
  
  #Set the new state for the next iteration
  x<-xnew
  
  #Save the state 
  xvec<-cbind(xvec,xnew)
}

#Prepare data for plotting
xvec<-t(xvec)
xvec<-as.data.frame(cbind(1:(layer_number-1),xvec))
colnames(xvec) <- c("time","heater 1","heater 2","heater 3","heater 4","heater 5","heater 6")

#Plot the trajectories using ggplot
xvec[-1,c("time","heater 1","heater 2","heater 3","heater 4","heater 5","heater 6")] %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=1)+
  theme_classic()+
  labs(x='time',y='heater temperature evolution')
