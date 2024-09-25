########################## ASL Week 3 ################################


cat("\014") 
graphics.off()  # clear all graphs
rm(list = ls()) # remove all files from your workspace

pacman::p_load(pacman, ggplot2 ) #installing neccesary packages 

########################## Exersice 0  ################################




# Parameters
n <- 360
p <- 1/18
lambda <- 20

# Define a range of values for plotting
x_binomial <- 0:40 # Binomial random variable (0 to  40)
x_poisson <- 0:40   # Poisson random variable (0 to 40)

# Compute the PDF for Binomial and Poisson distributions
pdf_binomial <- dbinom(x_binomial, size = n, prob = p)
pdf_poisson <- dpois(x_poisson, lambda = lambda)

# Create a plot
plot(x_binomial, pdf_binomial, type = "h", lwd = 2, col = "blue", ylim = c(0, max(pdf_binomial, pdf_poisson)), 
     xlab = "Value", ylab = "Probability", main = "Comparison of Binomial and Poisson Distributions")

# Add Poisson distribution to the plot
lines(x_poisson, pdf_poisson, type = "h", lwd = 2, col = "red")

# Add a legend
legend("topright", legend = c("Binomial (n=360, p=1/18)", "Poisson (λ=20)"), col = c("blue", "red"), lwd = 2)


### Takeaway: Poisson and Binomianl are alike when n is hig


########################## Exersice 1 ################################



## 

# Code provided by lecturer
##-----------------------------------------------------
## Earth quakes data
##-----------------------------------------------------


## Number of major world-wide earth quakes between 1900 and 2006 
ObsSeq <- c(13,14,8,10,16,26,32,27,18,32,36,24,22,23,22,18,25,21,
            21,14,8,11,14,23,18,17,19,20,22,19,13,26,13,14,22,24,
            21,22,26,21,23,24,27,41,31,27,35,26,28,36,39,21,17,22,
            17,19,15,34,10,15,22,18,15,20,15,22,19,16,30,27,29,23,
            20,16,21,21,25,16,18,15,18,14,10,15,8,15,6,11,8,7,18,
            16,13,12,13,20,15,16,12,18,15,16,13,15,16,11,11)
##-------------------------------------------------------
## EM algorithm for two component Poisson mixture model
##------------------------------------------------------
## Starting values
alphaEst <- 0.6
lamEst <- 2
etaEst <- 4

## Run EM algorithm
old_LL <-0
threshold <- 0.1
iter <- 0

##The function runs until the Log likelihood function stops changing. 
while (threshold > 1e-6){ 
  iter=iter+1
  xPrb <- alphaEst*dpois(ObsSeq,lamEst)/
    (alphaEst*dpois(ObsSeq,lamEst)+
       (1-alphaEst)*dpois(ObsSeq,etaEst))
  
  
  alphaEst <- sum(xPrb)/length(ObsSeq)
  
  lamEst <- sum(xPrb*ObsSeq)/sum(xPrb)
  
  etaEst <- sum((1-xPrb)*ObsSeq)/sum(1-xPrb)
  
  logLk <- sum(log(alphaEst*dpois(ObsSeq,lamEst)+
                     (1-alphaEst)*dpois(ObsSeq,etaEst)))
  cat("Iteration:",iter,
      "alpha:",alphaEst,"lam:",lamEst,"eta:",etaEst,
      "logLk",logLk,"\n")
  threshold <- abs(logLk-old_LL)
  old_LL <- logLk
}
##------------------------------------------------------
plot(ObsSeq,type="l")
points(1:length(xPrb),5+10*(1-xPrb),col="red",type="l")
##------------------------------------------------------
mx <- max(ObsSeq)
plot(1:mx,tabulate(ObsSeq),type="l",lwd=2)
points(1:mx,length(ObsSeq)*alphaEst*dpois(1:mx,lam=lamEst),
       col="red",type="l",lwd=2)
points(1:mx,length(ObsSeq)*(1-alphaEst)*dpois(1:mx,lam=etaEst),col="blue",
       type="l",lwd=2)
points(1:mx,
       length(ObsSeq)*(alphaEst*dpois(1:mx,lam=lamEst)+
                         (1-alphaEst)*dpois(1:mx,lam=etaEst)),
       col="green",type="l",lwd=2)



### c - reconcile with HMM book

#The graphs are the same





### d - reproduce the values for m=1 in table 1.2 and reproduce upper left plot in figure 1.4


#for m=1 just take the mean. 





########################## Exersice 2 ################################

#Loading data

dat <- read.table("C:/Users/au648278/OneDrive - Aarhus universitet/Dokumenter/PhD/Courses/ASL/Code/Week 3/YMisDat.txt")
## Construct Y array
Yarr <- array(0,dim=c(3,2,6))
for (rw in 1:36){
  i <- dat[rw,1] 
  j <- dat[rw,2]
  k <- dat[rw,3]
  Yarr[i,j,k] <- dat[rw,4]
}
## Construct W array
Warr <- array(1,dim=c(3,2,6))
Warr[1,1,1] <- 0
Warr[3,2,6] <- 0



## Lange Approach - example 9.4

two_way_ANOVA <- function(Y) {
  dims <- dim(Y)
  i <- dims[1]
  j <- dims[2]
  k <- dims[3]
  
  S <- apply(Y, 3, sum)
  
  alpha <- rowSums(S)  # Sum across rows, returns a vector of length i
  beta <- colSums(S)   # Sum across columns, returns a vector of length j
  
  mu <- sum(alpha) / (i * j * k)
  
  alpha <- (alpha / (j * k)) - mu
  beta <- (beta / (i * k)) - mu
  
  return(list(mu = mu, alpha = alpha, beta = beta))
}


MM_ANOVA <- function(Y, W) {
  dims <- dim(Y)
  i <- dims[1]
  j <- dims[2]
  k <- dims[3]
  
  mu <- 0
  alpha <- rep(0, i)
  beta <- rep(0, j)
  
  P <- array(0, dim = dim(Y)) # predicted values
  X <- array(0, dim = dim(Y))
  
  
  S=sum(Y,)
  
  old_rss <- Inf
  
  for (n in 1:100) { # mm iteration loop
    rss <- 0 # residual sum of squares
    
    for (ii in 1:i) {
      for (jj in 1:j) {
        for (kk in 1:k) {
          P[ii, jj, kk] <- mu + alpha[ii] + beta[jj]
          residual <- Y[ii, jj, kk] - P[ii, jj, kk]
          rss <- rss + W[ii, jj, kk] * residual^2
        }
      }
    }
    
    X <- W * Y + (1 - W) * P
    
    result <- two_way_ANOVA(X) # MM update
    mu <- result$mu
    alpha <- result$alpha
    beta <- result$beta
    
    if (abs(old_rss - rss) < 1e-8) { # check for convergence
      return(list(mu = mu, alpha = alpha, beta = beta))
    }
    
    old_rss <- rss
  }
  
  return(list(mu = mu, alpha = alpha, beta = beta))
}


result <- MM_ANOVA(Yarr, Warr)


