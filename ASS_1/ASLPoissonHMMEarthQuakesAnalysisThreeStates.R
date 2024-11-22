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
  #cat("Updated transition probability matrix:","\n")
  #print(TrnsPrb)
}


