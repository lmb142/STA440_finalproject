## metadata
#rm(list = ls())
#install.packages("devtools")
install.packages("coda")
library(coda)

## Standalone function
## Monte Carlo sampling event handler function
MC <- function(y, a = NULL, b = NULL, r = NULL, smc, prior, sampling, predictive = FALSE, summary = FALSE) {
  if (prior=="gamma" & sampling=="poisson" & predictive == FALSE){
    return (gammaPoisPosterior(y, a, b, smc))
  } else if (prior=="gamma" & sampling=="poisson" & predictive == TRUE){
    return (gammaPoisPredictive(y, a, b, smc))
  } else if (prior=="beta" & sampling=="bernoulli" & predictive == FALSE){
    return (betaBernoulliPosterior(y, a, b, smc))
  } else if (prior=="uniform" & sampling=="bernoulli" & predictive == FALSE){
  return (uniformBernoulliPosterior(y,smc))
  } else if (prior=="beta" & sampling=="negativebinomial" & predictive == FALSE){
    return (betaNegativeBinomialPosterior(y, a, b, r, smc))
  } else if (prior=="beta" & sampling=="geometric"  & predictive == FALSE){
    return(betaGeometricPosterior(y, a, b, smc))
  } else if (prior=="gamma" & sampling=="exponential" & predictive==FALSE){
    return(gammaExponentialPosterior(y, a, b, smc))
  } else{
    return("**Chosen distribution combinations are not supported**")
  }
}


## Normal prior on mean & Normal sampling distribution with variance known, mean posterior outputs
univariateNormMC_varKnown <- function(y, mu, tau, sigma, smc, output = NULL, summary = FALSE){
  varianceValue = 1 / (1/tau^2 + length(y) / sigma^2)
  meanValue = (mu /tau^2 + sum(y) / sigma^2) * varianceValue
  thetaValue <- rnorm(smc, meanValue, sqrt(varianceValue))
  predValue <- rnorm(smc, meanValue, sqrt(varianceValue + sigma^2))
  if (output=="mean"){
    return(thetaValue)
  } else if (output=="predictive"){
    return(predValue)
  } else{
    return(list(thetaValue, predValue))
  }
}


## Standalone function
## Normal prior on mean, Normal sampling distribution & inverse gamma prior on variance
## posterior outputs for mean & variance
univariateNormMC_varUnknown <- function(y, v0, sigma0, mu0, kappa0, smc, output = NULL, summary = FALSE){
  n <- length(y)
  vn <- v0 + n
  sigman <- 0.5 * ((length(y) - 1)*var(y) + v0*sigma0^2 + (n*kappa0)/(n + kappa0)*(mean(y)-mu0)^2)
  sigmaValue <- 1 / rgamma(smc, vn, sigman)
  thetaValue <- rnorm(smc, (kappa0*mu0 + n*mean(y))/(kappa0 + n), sqrt(sigmaValue/(kappa0 + n)))

  if (output=="mean"){
    return(thetaValue)
  } else if (output=="variance"){
    return(sigmaValue)
  } else{
    return(list(thetaValue, sigmaValue))
  }
}

## Standalone Function
## Gibbs sampler for Mean and Variance of Univariate Normal Distribution
univariateNorm_Gibbs <- function(y, v0, sigma0, mu0, tau0, smc, output=NULL, summary = FALSE){
  vn <- v0 + length(y)
  theta_0 <- mean(y)
  sigma_0 <- 1 / var(y)
  PHI <- matrix(nrow = smc, ncol=2)
  PHI[1,] <- phi <- c(theta_0, sigma_0)

  for (i in 2:smc){
    sigma_n <- (v0*sigma0^2 + sum(y-phi[1])^2)
    phi[2] <- rnorm(1, vn/2, sqrt(sigma_n/2))

    tau_n <- 1 / (1/tau^2 + length(y) / phi[2]^2)
    mu_n <- (mu0/tau0^2 + sum(y)/phi[2]) * tau_n
    phi[1] <- rnorm(1, mu_n, sqrt(tau_n))

    PHI[s,] <- phi
  }
  if (output=="both"){
    return (PHI)
  } else if (output=="mean"){
    return (PHI[,1])
  } else if (output=="variance"){
    return (PHI[,2])
  } else {
    return (PHI)
  }
}




## Standalone Function
## Note: this function uses logic from Peter Hoff's textbook: "A First Course in Bayesian
## Statistical Methods"
mvn_Gibbs <- function(y, mu0, Lambda0, v0, S0, smc, output=NULL, summary = FALSE) {
  n <- dim(Y)[1]
  y <- as.matrix(y)
  ybar <- apply(y, 2, mean)
  Sigma <- Cov(y)
  THETA <- SIGMA <-  NULL #PREDICTIVE

  for (i in 1:smc){
    Lambdan <- (solve(lambda0) + n * solve(Sigma))
    mun <- Lambdan%*% (solve(Lambda0)%*%mu0 + n*solve(Sigma)%*%ybar)
    theta <- rmvnorm(1, mun, Lambdan)

    S_theta <- (t(y)-c(theta))%*%t(t(y)-c(theta))
    Sigma <- solve(rwish, v0+n, solve(S0 + S_theta))

    predictive <- rmvnorm(1, theta, Sigma)

    THETA <- rbind(THETA, theta)
   # PREDICTIVE <- rbind(PREDICTIVE, predictive)
    SIGMA <- rbind(SIGMA, c(Sigma))
  }

  if (output=="mean"){
    return(THETA)
  } else if (output=="variance"){
    return(SIGMA)
 # } #else if (output=="predictive"){
    #return(PREDICTIVE)
  } else {
    return(list(THETA, SIGMA)) #PREDICTIVE
  }
}

simulationSummary <- function(values, type){
  if (type=="MC"){
    print(paste("Approximate expected value is: ", mean(x), sep=""))
    print(paste("Approximate variance value is: ", var(x), sep=""))
    print(paste("Monte Carlo standard error is: ", sd(x), sep=""))
    print(hist(x, main = "Approxmate Distribution", xlab = "Value", ylab = "Density"))
  } else if (type=="Gibbs") {
    print(paste("Approximate expected value is: ", mean(x), sep=""))
    print(paste("Approximate variance value is: ", var(x), sep=""))
    print(hist(x, main = "Approxmate Distribution", xlab = "Value", ylab = "Density"))
    print(paste("Autocorrelation is: ", acf(x), sep=""))
    print(paste("Effectivs sample size is: ", effectiveSize(x), sep=""))
  }
}

## The below functions are helper functions ------------------------------------

## Gamma prior & Pois sampling distribution with posterior outputs
gammaPoisPosterior <- function(y, a, b, smc){
  ret <- rgamma(smc, a + sum(y), b + length(y))
  return (ret)
}

## Gamma prior & Pois sampling distribution with predictive outputs
gammaPoisPredictive <- function(y, a, b, smc){
  theta <- rgamma(smc, a + sum(y), b + length(y))
  ret <- rpois(smc, theta)
  return (ret)
}

## Beta prior & Bernoulli sampling distribution with posterior outputs
betaBernoulliPosterior <- function(y, a, b, smc) {
  ret <- rbeta(smc, a + sum(y), b + length(y) - sum(y))
  return (ret)
}

## Unif prior & Bernoulli sampling distribution with posterior outputs
uniformBernoulliPosterior <- function(y, smc) {
  ret <- rbeta(smc, sum(y), length(y) - sum(y))
  return (ret)
}

## Beta prior & Negative Binomial sampling distribution with posterior outputs
betaNegativeBinomialPosterior <- function(y, a, b, r, smc){
  ret <- rbeta(smc, a + r*length(y), b + sum(y))
  return (ret)
}

## Beta prior & Gemeotric sampling distribution with posterior outputs
betaGeometricPosterior <- function(y, a, b, smc){
  ret <- rbeta(smc, a + n, b + sum(y))
  return (ret)
}

## Gamma prior & Exponential sampling distribution with posterior outputs
gammaExponentialPosterior <- function(y, a, b, smc){
  ret <- rgamma(smc, a + length(y), b + sum(y))
  return (ret)
}

## Note: this helper function is taken from Peter Hoff's textbook: "A First Course in Bayesian
## Statistical Methods"
rwish<-function(n,nu0,S0){
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n){
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}


