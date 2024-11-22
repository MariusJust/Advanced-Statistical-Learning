# ADVANCED STATISTICAL LEARNING
# Name: HMMPoisExpectationsFct.R
# Author: Asger Hobolth
# Purpose:
# Calculates expected transition count matrix
# Calculates likelihood
# Calculates posterior probability
#--------------------------------------------------------
# Input:
# InitProb: Initial probabilities
# TransProb: Transition probabilities;
#            Probabilities between hidden states
#            nHS times nHS matrix
# Lambda: Emission probabilities for Poisson 
#         Vector of length nHS
# ObsSeq: Observed sequence
#
# Output:
# $TransCnt: Expected transition count matrix
# $PostProb: Posterior probability matrix
# $logLk: Log likelihood
#---------------------------------------------------------
HMMPoisExpectationsFct <- function(InitProb,TransProb,Lambda,ObsSeq){
  nObs <- length(ObsSeq)
  nHS <- nrow(TransProb)
  #-------------------
  # Forward algorithm 
  #-------------------
  # Define Forward matrix
  ForwardMat <- matrix(0,nrow=nObs,ncol=nHS)
  # Start conditions
  ForwardMat[1,] <- InitProb*dpois(ObsSeq[1],Lambda)
  # Determine ForwardMat by recursion
  for (k in 2:nObs){
    for (j in 1:nHS){
      ForwardMat[k,j] <- dpois(ObsSeq[k],Lambda[j])*
                         sum(TransProb[,j]*
                             ForwardMat[k-1,])
    }
  }
  ForwardLik <- sum(ForwardMat[nObs,])
  #cat("log-Likelihood from Forward algorithm:",log(ForwardLik),"\n")
  #--------------------
  # Backward algorithm
  #--------------------
  # Define BackwardMat
  BackwardMat <- matrix(0,nrow=nObs,ncol=nHS)
  # Start condition
  BackwardMat[nObs,] <- rep(1,nHS)
  # Determine BackwardMat by recursion
  for (k in (nObs-1):1){
    for (j in 1:nHS){
      BackwardMat[k,j] <- sum(TransProb[j,1:nHS]*
                              dpois(ObsSeq[k+1],Lambda)*  
                              BackwardMat[k+1,])
    }
  }
  BackwardLik <- sum(InitProb*
                        dpois(ObsSeq[1],Lambda)*
                        BackwardMat[1,])
  #cat("log-Likelihood from Backward algorithm:",log(BackwardLik),"\n")
  ##-----------------------
  ## Posterior probability
  ##-----------------------
  PostProb <- exp(log(BackwardMat)+log(ForwardMat)-log(BackwardLik))
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
  }
  output <- list()
  output$TransCnt <- TransCnt
  output$PostProb <- PostProb
  output$logLk <- log(BackwardLik)
  return(output)
}

