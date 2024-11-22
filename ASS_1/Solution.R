#################################### Mandatory assignment ########################################


## cleanup 

cat("\014") 
graphics.off()  # clear all graphs
rm(list = ls()) # remove all files from your workspace
pacman::p_load(pacman, ggplot2, dplyr ) #installing neccesary packages 


################################### Part 1 #######################################################
#                                                                                                #                
#                       Parameter estimation using EM                                            #
#                                                                                                #
##################################################################################################

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
## Plot data
plot(ObsSeq,type="l")

HMMPoisExpectationsFct <- function(InitProb,TransProb,Lambda,ObsSeq){
  
  nObs <- length(ObsSeq) # number of observed states
  nHS <- nrow(TransProb) # Number of hidden states 
  
  
  #-------------------
  # Forward algorithm 
  #-------------------
  
  # Define Forward matrix
  ForwardMat <- matrix(0,nrow=nObs,ncol=nHS)
  # Start conditions
  ForwardMat[1,] <- InitProb*dpois(ObsSeq[1],Lambda) # initialising step as in proposition 2 (page 66) in the book
  # Determine ForwardMat by recursion
  #Itterating over the number of observations and number of hidden states. At each iteration the recursion (equation 4.3 in the book) is calculated
  for (k in 2:nObs){
    for (j in 1:nHS){
      ForwardMat[k,j] <- dpois(ObsSeq[k],Lambda[j])* # probability of observing variable k, given a poisson distribution with frequency lambda. 
        sum(TransProb[,j]* 
              ForwardMat[k-1,]) #summing over all the transformation probabilities (going from state i to j - for all i), times the previous valuesof the forward matrix
    }
  }
  ForwardLik <- sum(ForwardMat[nObs,]) # The forward likelihood is calculated as the sum of the last column in the forward matrix
  #cat("log-Likelihood from Forward algorithm:",log(ForwardLik),"\n")
  #--------------------
  # Backward algorithm
  #--------------------
  # Define BackwardMat
  BackwardMat <- matrix(0,nrow=nObs,ncol=nHS)
  # Start condition
  BackwardMat[nObs,] <- rep(1,nHS) # all 3 hidden states are initialised with 1.  
  # Determine BackwardMat by recursion
  for (k in (nObs-1):1){
    for (j in 1:nHS){
      BackwardMat[k,j] <- sum(TransProb[j,1:nHS]*
                                dpois(ObsSeq[k+1],Lambda)*  
                                BackwardMat[k+1,])
    } # at each backward step we take the sum over the transformation probability times the probability of observing the k+1'th variable, given it is poisson distributed. 
  } # Multiplied with the previous input of the backward matrix 
  BackwardLik <- sum(InitProb*
                       dpois(ObsSeq[1],Lambda)*
                       BackwardMat[1,]) #The likelihood is the sum of the initial probabiliy times the probability of drawing the first observation from a poisson distribution, times the first column of the backward matrix 
  #cat("log-Likelihood from Backward algorithm:",log(BackwardLik),"\n")
  ##-----------------------
  ## Posterior probability
  ##-----------------------
  PostProb <- exp(log(BackwardMat)+log(ForwardMat)-log(BackwardLik)) # avoid underflow issues by exponentiating - formula is equivalevnt to proposition 5 (4.10)
  ##----------------------------
  ## Expected transition counts
  ##----------------------------
  TransCnt <- matrix(0,nrow=nHS,ncol=nHS)
  for (i in 1:nHS){
    for (j in 1:nHS){
      
      
      Probij <- ForwardMat[1:(nObs-1),i]*BackwardMat[2:nObs,j]/
        BackwardLik*TransProb[i,j]*dpois(ObsSeq[2:nObs],Lambda[j])
      
      TransCnt[i,j] <- sum(Probij)
    }
  } # calculating the expected transistion counts. Similar to formula 4.11 in the book 
  output <- list()
  output$TransCnt <- TransCnt
  output$PostProb <- PostProb
  output$logLk <- log(BackwardLik)
  return(output)
} 




##--------------------------------------------------
## Initial values for the EM with three hidden states
##--------------------------------------------------
## Three hidden states
## Initial probabilities
IntPrb <- rep(1,3)/3
## Transition probabilities
TrnsPrb <- matrix(c(0.8,0.1,0.1,
                    0.1,0.8,0.1,
                    0.1,0.1,0.8),
                  byrow=TRUE,nrow=3,ncol=3)
## Rates for emission
Lam <- c(10,20,30)
## Iteration 0 in Table 4.2:
HMMexpct <- HMMPoisExpectationsFct(IntPrb,TrnsPrb,Lam,ObsSeq)
cat("log-Likelihood:",HMMexpct$logLk,"\n")


##------------------------------------------------
## Estimate transition matrix and rates using EM
##-----------------------------------------------
## Number of iterations
nIter <- 20
## EM algorithm
for (iter in 1:nIter){
  #cat("--------------------------------------","\n")
  HMMexpct <- HMMPoisExpectationsFct(IntPrb,TrnsPrb,Lam,ObsSeq)
  #cat("\n")
  cat("Iteration:",iter,"log-likelihood:",HMMexpct$logLk,"\n")
  ## Updating initial probability
  IntPrb <- HMMexpct$PostProb[1,]
  ## Updating transition probability matrix
  TrnsPrb <- HMMexpct$TransCnt/rowSums(HMMexpct$TransCnt)
  ## Updating rates in Poisson distribution 
  for (i in 1:nrow(TrnsPrb)){
    Lam[i] <- sum( HMMexpct$PostProb[,i]*ObsSeq )/sum( HMMexpct$PostProb[,i] )
  }
  cat("Updated transition probability matrix:","\n")
  print(TrnsPrb)

}







################################### Part 2 #######################################################
#                                                                                                #                
#                             Posterior states and decoding                                      #
#                                                                                                #
##################################################################################################

#Fetching the posterior probabilities from the already existing function
post_prob_state_1 <- HMMexpct$PostProb[, 1]
post_prob_state_2 <- HMMexpct$PostProb[, 2]
post_prob_state_3 <- HMMexpct$PostProb[, 3]


# Set up a 3-row plotting space
par(mfrow=c(3,1))

# Plot for State 3
plot(post_prob_state_3, type="h", lwd=2, ylim=c(0,1), 
     main="State 3", xlab="Time", ylab="Posterior Probability")

# Plot for State 2
plot(post_prob_state_2, type="h", lwd=2, ylim=c(0,1), 
     main="State 2", xlab="Time", ylab="Posterior Probability")

# Plot for State 1
plot(post_prob_state_1, type="h", lwd=2, ylim=c(0,1), 
     main="State 1", xlab="Time", ylab="Posterior Probability")



########### Decoding 


####### Posterior ( local decoding )


# I now have all the posterior probability, and therefore I simply choose the state that has the highest posterior probability


posterior_decoding <- function(PostProb) {
  # For each time step, select the state with the highest posterior probability
  decoded_states <- apply(PostProb, 1, which.max)
  return(decoded_states)  # This gives the most probable state at each time step
}

decoded_states <- posterior_decoding(HMMexpct$PostProb)

# Observed sequence of data
plot(ObsSeq, type = "l", xlab = "Year", ylab = "Earthquakes", main = "Figure 5.5", ylim = c(0, 50))

# Add horizontal lines for each state's mean (lambda values)
abline(h = Lam[1], col = "red", lty = 2)  # Mean of state 1
abline(h = Lam[2], col = "blue", lty = 2) # Mean of state 2
abline(h = Lam[3], col = "green", lty = 2) # Mean of state 3

# Highlight periods where the decoded state is active
for (i in 1:length(ObsSeq)) {
  if (decoded_states[i] == 1) {
    points(i, ObsSeq[i], col = "red", pch = 19)  # Highlight state 1
  } else if (decoded_states[i] == 2) {
    points(i, ObsSeq[i], col = "blue", pch = 19) # Highlight state 2
  } else if (decoded_states[i] == 3) {
    points(i, ObsSeq[i], col = "green", pch = 19) # Highlight state 3
  }
}


######## Global (Viterbi)


# Viterbi decoding function using Poisson HMM
pois.HMM.viterbi <- function(x, mod) {
  n <- length(x)
  xi <- matrix(0, n, mod$m)
  
  # Initial probabilities (delta * Poisson likelihood)
  foo <- mod$delta * dpois(x[1], mod$lambda)
  xi[1, ] <- foo / sum(foo)
  
  # Recursive Viterbi algorithm
  for (i in 2:n) {
    foo <- apply(xi[i - 1, ] * mod$gamma, 2, max) * dpois(x[i], mod$lambda)
    xi[i, ] <- foo / sum(foo)
  }
  
  # Backtrack to find the most probable states
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for (i in (n - 1):1) {
    iv[i] <- which.max(mod$gamma[, iv[i + 1]] * xi[i, ])
  }
  return(iv)
}


# Define the model structure to pass to Viterbi decoding
mod <- list(m = 3,                # Number of hidden states
            delta = IntPrb,        # Initial probabilities
            gamma = TrnsPrb,       # Transition matrix
            lambda = Lam)          # State-dependent Poisson means


# Perform Viterbi decoding
viterbi_states <- pois.HMM.viterbi(ObsSeq, mod)



# Plot the observed data
plot(ObsSeq, type = "l", xlab = "Year", ylab = "Earthquakes", main = "Figure 5.6", ylim = c(0, 50))

# Add horizontal lines for each state's mean (lambda values)
abline(h = Lam[1], col = "red", lty = 2)  # Mean of state 1
abline(h = Lam[2], col = "blue", lty = 2) # Mean of state 2
abline(h = Lam[3], col = "green", lty = 2) # Mean of state 3

# Highlight periods wwe here the Viterbi-decoded state is active
for (i in 1:length(ObsSeq)) {
  if (viterbi_states[i] == 1) {
    points(i, ObsSeq[i], col = "red", pch = 19)  # Highlight state 1
  } else if (viterbi_states[i] == 2) {
    points(i, ObsSeq[i], col = "blue", pch = 19) # Highlight state 2
  } else if (viterbi_states[i] == 3) {
    points(i, ObsSeq[i], col = "green", pch = 19) # Highlight state 3
  }
}









################################### Part 3 #######################################################
#                                                                                                #                
#                               State prediction                                                 #
#                                                                                                #
##################################################################################################



################################### Part 4 #######################################################
#                                                                                                #                
#                                 Forecasting                                                    #
#                                                                                                #
##################################################################################################


################################### Part 5 #######################################################
#                                                                                                #                
#                                Model Selection                                                 #
#                                                                                                #
##################################################################################################


################################### Part 6 #######################################################
#                                                                                                #                
#                                 Confidence Intervals                                           #
#                                                                                                #
##################################################################################################


# Function to estimate parameters using EM algorithm
estimate_HMM_params <- function(init_probs, trans_probs, lambdas, obs_seq, n_iter = 20) {
  for (iter in 1:n_iter) {
    HMMexpct <- HMMPoisExpectationsFct(init_probs, trans_probs, lambdas, obs_seq)
    init_probs <- HMMexpct$PostProb[1,]
    trans_probs <- HMMexpct$TransCnt / rowSums(HMMexpct$TransCnt)
    for (i in 1:nrow(trans_probs)) {
      lambdas[i] <- sum(HMMexpct$PostProb[, i] * obs_seq) / sum(HMMexpct$PostProb[, i])
    }
  }
  return(list(init_probs = init_probs, trans_probs = trans_probs, lambdas = lambdas))
}



# Function from the book to generate bootstrap sample
pois.HMM.generate_sample <- function(ns, mod) {
  mvect <- 1:3  # m is the number of states
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob = mod$delta)
  
  # Generate the state sequence
  for (i in 2:ns) {
    state[i] <- sample(mvect, 1, prob = mod$gamma[state[i - 1], ])
  }
  
  # Generate the observation sequence based on the state
  x <- rpois(ns, lambda = mod$lambda[state])
  
  return(x)
}




# Fit the model to get MLE estimates (initial fit)
MLE_fit <- estimate_HMM_params(IntPrb, TrnsPrb, Lam, ObsSeq)



# Number of bootstrap iterations
B <- 1000  # Adjust this value as needed

# Store the bootstrap parameter estimates
bootstrap_estimates <- matrix(NA, nrow = B, ncol = 15)  # Assuming 3 lambdas, 9 gamma elements, 3 deltas

# Bootstrap process
set.seed(123)  # For reproducibility
for (b in 1:B) {
  # Step 2a: Generate a bootstrap sample
  bootstrap_sample <- pois.HMM.generate_sample(length(ObsSeq), MLE_fit)
  
  # Step 2b: Re-estimate parameters using EM algorithm on the bootstrap sample
  bootstrap_fit <- estimate_HMM_params(MLE_fit$init_probs, MLE_fit$trans_probs, MLE_fit$lambda, bootstrap_sample)
  
  # Step 2c: Record the parameter estimates
  bootstrap_estimates[b, 1:3] <- bootstrap_fit$lambdas
  bootstrap_estimates[b, 4:12] <- as.vector(bootstrap_fit$trans_probs)
  bootstrap_estimates[b, 13:15] <- as.vector(bootstrap_fit$init_probs)
}

# Column names for the parameters (3 lambdas, 9 transition probabilities)
colnames(bootstrap_estimates) <- c("lambda1", "lambda2", "lambda3", 
                                   "gamma11", "gamma12", "gamma13", 
                                   "gamma21", "gamma22", "gamma23", 
                                   "gamma31", "gamma32", "gamma33",
                                   "delta1", "delta2", "delta3")

# Compute the 90% confidence intervals
conf_intervals <- apply(bootstrap_estimates, 2, function(x) quantile(x, probs = c(0.1, 0.9)))

# Convert to a data frame for display
conf_intervals_df <- as.data.frame(t(conf_intervals))
colnames(conf_intervals_df) <- c("10th Percentile", "90th Percentile")

# Add the MLE estimates to the table
conf_intervals_df$MLE <- colMeans(bootstrap_estimates)

# Print the results
# Prevent scientific notation and print results with 1 decimal place
options(scipen = 999)

# Round to 1 decimal place before printing

conf_intervals_df %>% mutate_if(is.numeric, ~ scales::number(., accuracy = 0.001))
# Print the data frame
print(conf_intervals_df)
















