library(Rlab)

#### Calculating first and second moments of a vector ####

## Compute first moment for the sample around origin -> Sample_Mean
first_moment <- function(vec){
  n <- length(vec)
  sum_total <- sum(vec, na.rm=TRUE)
  return (sum_total/n)
}

## Compute second moment for the sample around the sample mean -> Sample_Variance ##
second_moment <- function(vec){
  mu <- first_moment(vec)
  n <- length(vec)
  square_total <- sum((vec - mu)^2, na.rm=TRUE)
  return (square_total/n)
}

#################################################
#### Method of Moment for each distribution ####
###############################################
poisson <- function(vec){
  lambda_hat = first_moment(vec)
  return(lambda_hat)
}

gamma <- function(vec){
  
  mu = first_moment(vec)
  var = second_moment(vec)
  alpha_hat <- (mu**2)/var
  beta_hat <- var/mu
  list(alpha_hat = alpha_hat, beta_hat = beta_hat)
}

beta <- function(vec){
  mu = first_moment(vec)
  var = second_moment(vec)
  temp = (mu*(1-mu)/var) - 1
  
  alpha_hat = mu * temp
  beta_hat = (1-mu) * temp
  
  list(alpha_hat = alpha_hat, beta_hat = beta_hat)
}

geometric <- function(vec){
  p_hat = 1/first_moment(vec) -1
  return(p_hat)
}

uniform <- function(vec){
  mu =  first_moment(vec)
  var = second_moment(vec)
  a_hat = mu - (var*sqrt(3))
  b_hat = mu + (var*sqrt(3))
  list(a_hat = a_hat, b_hat = b_hat)
}

exponential <- function(vec){
  beta_hat = 1/first_moment(vec)
  return(beta_hat)
}



chi_squared <- function(vec){
  p_hat = first_moment(vec)
  return(p_hat)
}

t <- function(vec){
  mu = first_moment(vec)
  if(mu == 0) 
    return(0)
  else{
    v_hat = 2*mu/(mu-1)
    return(v_hat)
  }
    
}

normal <- function(vec){
  
  mu_hat <- first_moment(vec)
  sd_hat <- sqrt(second_moment(vec))
  
 list(mu_hat = mu_hat, sd_hat = sd_hat)
  
}

binomial <- function(vec){
  
  mu <- first_moment(vec)
  var <- second_moment(vec)
  
  n_hat <- mu/(1-(var/mu))
  p_hat <- 1 - (var/mu)
  
  list(n_hat = n_hat, p_hat = p_hat)
}

bernoulli <- function(vec){
  
  p_hat = first_moment(vec)
  return(p_hat)
}

###############################################################################
#### Function to generate data of distribution and estimate the parameters ####
###############################################################################

method_of_moment <- function(n, distribution){
  
  cat(paste("Distribution :", distribution ,"\n" ))
  if( distribution == "poisson" ){
    lambda = 0.01
    vec = rpois(n, lambda)
    lambda_hat = poisson(vec)
    print(paste0("Population Parameter:   ", "lambda = ", lambda))
    print(paste0("Estimated parameter:", "lambda_hat = ",round(lambda_hat,4) ))
  }
  
  if( distribution == "gamma" ){
    alpha = 2
    beta = 0.4
    vec <- rgamma(n, alpha, scale = beta) 
    l <- gamma(vec)
    print(paste0("Population Parameters:    ", "alpha = ", alpha, "               beta = ", beta))
    print(paste0("Estimated parameters: ", "alpha_hat = ",round(l$alpha_hat,4), "      beta_hat = ", round(l$beta_hat,4) ) )
  }
  
  if( distribution == "beta" ){
    alpha = 3
    beta = 5
    vec <- rbeta(n, alpha, beta) 
    l <- beta(vec)
    print(paste0("Population Parameters:    ", "alpha = ", alpha, "               beta = ", beta))
    print(paste0("Estimated parameters: ", "alpha_hat = ",round(l$alpha_hat,4), "      beta_hat = ", round(l$beta_hat,4) ) )
  }
  
  if( distribution == "normal" ){
    mu = 10
    sd = 9
    vec <- rnorm(n, mu, sd) 
    l <- normal(vec)
    print(paste0("Population Parameters:    ", "mu = ", mu, "               sd = ", sd))
    print(paste0("Estimated parameters: ", "mu_hat = ",round(l$mu_hat,4), "      sd_hat = ", round(l$sd_hat,4) ) )
  }
  
  if( distribution == "geometric" ){
    p = 0.7
    vec = rgeom(n, p)
    p_hat = geometric(vec)
    print(paste0("Population Parameter:   ", "p = ", p))
    print(paste0("Estimated parameter:", "p_hat = ",round(p_hat,4) ))
  }
  
  if( distribution == "uniform" ){
    a = 1
    b = 5
    vec <- runif(n, a, b) 
    l <- uniform(vec)
    print(paste0("Population Parameters:    ", "a = ", a, "               b = ", b))
    print(paste0("Estimated parameters: ", "a_hat = ",round(l$a_hat,4), "      b_hat = ", round(l$b_hat,4) ) )
  }
  
  if( distribution == "exponential" ){
    beta = 5
    vec = rexp(n, beta)
    beta_hat = exponential(vec)
    print(paste0("Population Parameter:   ", "beta = ", beta))
    print(paste0("Estimated parameter:", "beta_hat = ",round(beta_hat,4) ))
  }
  
  if( distribution == "chi_squared" ){
    p = 3
    vec = rchisq(n, p)
    p_hat = chi_squared(vec)
    print(paste0("Population Parameter:   ", "p = ", p))
    print(paste0("Estimated parameter:", "p_hat = ",round(p_hat,4) ))
  }
  
  
  if( distribution == "binomial" ){
    n_ = 10
    p = 0.7
    vec <- rbinom(n, n_, p) 
    l <- binomial(vec)
    print(paste0("Population Parameters:    ", "n_ = ", n_, "               p = ", p))
    print(paste0("Estimated parameters: ", "n_hat = ",round(l$n_hat,4) , "      p_hat = ", round(l$p_hat,4)) )
  }
  
  if( distribution == "bernoulli" ){
    p = 0.6
    vec = rbern(n, p)
    p_hat = bernoulli(vec)
    print(paste0("Population Parameter:   ", "p = ", p))
    print(paste0("Estimated parameter:", "p_hat = ", round(p_hat,4) ))
  }
  
  if( distribution == "t" ){
    v = 3
    vec = rt(n, df = v)
    v_hat = t(vec)
    print(paste0("Population Parameter:   ", "v = ", v))
    print(paste0("Estimated parameter:", "v_hat = ",round(v_hat,4) ))
  }
  
  cat("\n")
}


################################################################
#### Main wrapper function to run mom for all distributions ####
################################################################
main <- function(){
  
  distributions <- c("bernoulli", "binomial", "geometric", "poisson", "uniform", "normal", "exponential", "gamma", "beta", "t", "chi_squared")
  
  n = 100
  for (distr in distributions){
    method_of_moment(n, distr)
  }
  
}


#### call to the main function to run mom ####
main()
 
method_of_moment(100, "exponential")

