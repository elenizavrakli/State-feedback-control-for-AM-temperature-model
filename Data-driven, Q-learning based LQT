library(Matrix) 
library(e1071)
library(pracma)
library(devtools)
library(expm)
library(MASS)
library(tidyverse)

#Setting the seed so that the results are reproducible 
set.seed(10)

n=6 #number of states 
m=7 #number of inputs

#Choose number of data to be simulated, sample size and iterations for training
n_dat <- 5000
n_sampLS <-3000
n_iter <- 100

#Initialise dimensions for all vectors and matrices that will be used for data generation 
x <- matrix(rep(0,n*(n_dat+1)),nrow=n) #matrix storing all state vectors
u <- matrix(rep(0,m*(n_dat+1)),nrow=m) #matrix storing all input vectors
r <- matrix(rep(0,n*(n_dat+1)),nrow=n) #reference signal
xaug <- matrix(rep(0,(n+n)*(n_dat+1)),nrow=n+n) #augmented state including the reference
zaug <- matrix(rep(0,(n+n+m)*(n_dat+1)),nrow=n+n+m) #augmented state including reference and input

#Initialise dimensions for all vectors and matrices that will be used during training 
lhsvi <- matrix(rep(0,(m+n+n)^2*(n_dat)),nrow=n_dat) #Kronecker product vectors
uest <- array(rep(0,n_iter*m*(n_dat+1)),c(m,n_iter,n_dat+1)) #estimate of optimal control during training
xiaug <- matrix(rep(0,(n+n+m)*(n_dat+1)),nrow=n+n+m) #augmented state includign state, reference and control estimate
costvi <- matrix(rep(0,(n_dat+1)),nrow=n_dat+1) #estimated cost during training, using the control estimate for the future cost term
Kgain <- array(rep(0,m*(n+n)*(n_iter+1)),c(n_iter+1,m,n+n))  #estimates of optimal gain

mu <- 0.001 #refularisation term for Least Squares

errnorm <- matrix(rep(0,n_iter*(n_dat)),nrow=n_iter) #stopping rule for training

noisevar <- 0.5 #variance for the noise to the input to make it persistently exciting

#Define system matrices for data simulation
#State matrix
A <- matrix(c(0.992,0.0018,0,0,0,0,
              0.0023,0.9919,0.0043,0,0,0,
              0,-0.0042,1.0009,0.0024,0,0,
              0,0,0.0013,0.9979,0,0,
              0,0,0,0,0.9972,0,
              0,0,0,0,0,0.9953),
            nrow = n, byrow = TRUE)

#Input matrix
B <- matrix(c(1.0033,0,0,0,0,0,-0.2175,
              0,1.0460,0,0,0,0,-0.0788,
              0,0,1.0326,0,0,0,-0.0020,
              0,0,0,0.4798,0,0,-0.0669,
              0,0,0,0,0.8882,0,0.1273,
              0,0,0,0,0,1.1699,-0.1792),
            nrow = n, ncol = m, byrow = TRUE)

#Matrix generating the reference signal 
Ft <- diag(n)

#Augmented system state and input matrices
Tmatr <-  as.matrix(bdiag(A,Ft))
B1 <- as.matrix(rbind(B,matrix(rep(0,n*m),n,m)))

#Weighting matrices for performance index
Q <- diag(n)
C1 <- as.matrix(cbind(diag(n),-diag(n)))
Q1 <- t(C1)%*%Q%*%C1
R <- diag(m)

gam <- 0.9 #discount factor

#Gain used for data generation
Kgain[1,,]<- matrix(c(  0.7395,   -0.0076,   -0.0003,   -0.0264,    0.0194,   -0.0170,   0.00000001,0.0000000033,-0.00000000225,0.000000000934,-0.0000000137, -0.0000000008432,
                        -0.0076,    0.7430 ,   0.0031,   -0.0093,    0.0068,   -0.0060, -0.00000011, 0.000000621, 0.0000000341,-0.000000000163, 0.000000009215,-0.000000000132,
                        -0.0003,   -0.0033,    0.7599,    0.0021,    0.0002,   -0.0002,  0.00000000268,-0.0000001, 0.000000471, 0.00000000141,-0.0000000004333,0.0000000092,
                        -0.0126,   -0.0042,    0.0016,    1.0971,    0.0092,   -0.0079,  0.000000032, 0.00000000012683, -0.00000074, -0.0000000031,0.0000000024,0.00000000263,
                        0.0171,    0.0058,    0.0002,    0.0170,    0.8179,    0.0108,   -0.00000000227, 0.000000000953, -0.00000000232, 0.00000000009, -0.00000001, 0.000000001,
                        -0.0198,   -0.0067,   -0.0002,   -0.0193,    0.0143,    0.6823,  -0.000000001, 0.0000000421, 0.00000000002001, -0.0000004033, -0.0000000002321, 0.0000000123,
                        -0.1525,   -0.0519,   -0.0018,   -0.1412,    0.1091,   -0.0977,  0.00000083, 0.000000000563, 0.000000135,0.00000000131, -0.000000001,0.000000009921),nrow=m,byrow=TRUE)

#Initialise estimate of kernel matrix H 
Holdvi <- diag(n+n+m)

#Initialise values of vectors used during data generation 
x[,1]<- rep(50,n)
r[,1]<-c(190,180,170,165,160,155)
u[,1] <- sqrt(noisevar)*rnorm(m)
xaug[,1] <- c(x[,1],r[,1])
zaug[,1] <- c(xaug[,1],u[,1])

#-----------------Data Generation-----------------------#
for(k in 1:n_dat){
  x[,k+1] = A%*%x[,k] + B%*%u[,k]
  r[,k+1] = r[,k]
  
  u[,k+1] = -Kgain[1,,]%*%xaug[,k+1]+sqrt(noisevar)*rnorm(m)+10*sin(runif(m)*pi*k/5)+8*sin(runif(m)*2*pi*k/5)+7*sin(runif(m)*3*pi*k/5)+6*sin(runif(m)*4*pi*k/5)+4*sin(runif(m)*5*pi*k/5)+3*sin(runif(m)*6*pi*k/5)+0.5*sin(runif(m)*7*pi*k/5)
  
  xaug[,k+1] = c(x[,k+1],r[,k+1])
  zaug[,k+1] = c(xaug[,k+1],u[,k+1])
}

#-------------Kronecker Product terms-------------------#
for(k in 1:n_dat){
  lhsvi[k,] = kronecker(t(zaug[,k]),t(zaug[,k]))
}

#---------------------Training--------------------------#
for(i in 1:n_iter){
  uest[,i,1] =  -Kgain[i,,]%*%xaug[,1] #estimate optimal control using current estimate of the optimal gain
  xiaug[,1] = c(xaug[,1],uest[,i,1]) #augment state with control estimate, to be used in cost term
  
  for(k in 1:n_dat){
    uest[,i,k+1] = -Kgain[i,,]%*%xaug[,k+1] #estimate optimal control using current estimate of the optimal gain
    xiaug[,k+1] = c(xaug[,k+1],uest[,i,k+1]) #augment state with control estimate, to be used in cost term
    #cost term including current cost + expected future cost based on the optimal control estimate
    costvi[k] = t(xaug[,k])%*%Q1%*%xaug[,k]+t(u[,k])%*%R%*%u[,k]+gam*t(xiaug[,k+1])%*%Holdvi%*%xiaug[,k+1]
  }
  #Initialise and calculate matrix and vector for Least squares
  summat = matrix(rep(0,(n+n+m)^4),nrow=(n+n+m)^2)
  sumvec = matrix(rep(0,(n+n+m)^2),nrow=(n+n+m)^2)
  for(l in 1:(n_sampLS)){
    summat = summat + lhsvi[l,]%*%t(lhsvi[l,])
    sumvec = sumvec + lhsvi[l,]*costvi[l]
  }
  
  leftmat = summat + mu*diag((n+n+m)^2) #add regularisation term
  rightvec = sumvec
  
  HLS_est = ginv(leftmat)%*%rightvec #estimate of kernel matrix obtained form LS
  
  errnorm[i] = norm(matrix(HLS_est,nrow=(n+n+m))-Holdvi) #check whether estimation has converged
  if(errnorm[i]<=0.001){
    break
  }
  
  #Use kernel matrix estimate to obtain the control gain estimate
  Holdvi = matrix(HLS_est,nrow=n+n+m)
  Hxx = Holdvi[1:(n+n),1:(n+n)]
  Hxu = Holdvi[1:(n+n),(n+n+1):(n+n+m)]
  Hux = Holdvi[(n+n+1):(n+n+m),1:(n+n)]
  Huu = Holdvi[(n+n+1):(n+n+m),(n+n+1):(n+n+m)]
  Kgain[i+1,,] = ginv(Huu)%*%Hux
}

#-----Simulate system behaviour with estimated gain-----#

#Initialise
x <- rep(50,n)
xopt <- c(190,180,170,165,160,155)
X<-c(x,xopt)
Kest=Kgain[i,,] #set the control gain equal to the final estimate from training
layer_number = 100 #how long to simulate the behaviour for

xvec <- x #vector where all states are saved to obtain the trajectories

for(layr in 1:layer_number){
  #Calculate the optimal control
  u <- - Kest%*%X
  
  #Simulation from the state space system, applying the optimal control 
  xnew<-Tmatr%*%X+B1%*%u
  X<-xnew
  
  #Save the new state
  xvec<-cbind(xvec,xnew[1:n])
}

#Prepare data for plotting
xvec<-t(xvec)
xvec<-as.data.frame(cbind(1:(layer_number+1),xvec))
colnames(xvec) <- c("time","heater 1","heater 2","heater 3","heater 4","heater 5","heater 6")

#Plot the trajectories using ggplot
xvec[,c("time","heater 1","heater 2","heater 3","heater 4","heater 5","heater 6")] %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=1)+
  theme_classic()+
  labs(x='time',y='heater temperature evolution')

