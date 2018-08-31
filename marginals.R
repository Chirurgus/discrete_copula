# Created by Oleksandr Sorochynskyi
# On 12 06 18

margin_empir <- function(x) {
  F <- ecdf(x);
  list(param= NA, cdf= F, df= NULL, ran= NULL);
}

margin_poiss <- function(x) {
  lambda <- mean(x);
  cdf <- function(t) ppois(t, lambda);
  df <- function(t) dpois(t, lambda);
  ran <- function(t) rpois(t, lambda);
  list(param= lambda, cdf= cdf, df= df, ran= ran);
}

margin_nbinom <- function(x) {
  dl.dr <- function(r, vec) {
    nobs <- length(vec)
    sum(digamma(vec+r)) - nobs*digamma(r) + nobs*log(r/(r+sum(vec)/nobs))
  }
  
  #mle 
  r_hat <- uniroot(function(r) dl.dr(r,x), interval = c(0.00001, 1000),tol=0.0001)$root
  p_hat <- { s <- sum(x); s/(length(x)*r_hat + s) }
  
  cdf <- function(t) pnbinom(t,r_hat,p_hat);
  df <- function(t) dnbinom(t,r_hat,p_hat);
  ran <- function(t) rnbinom(t, r_hat,p_hat);
  list(param= c(r_hat,p_hat), cdf= cdf, df= df, ran= ran);
}

