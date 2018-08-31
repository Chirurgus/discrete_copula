# Created by Oleksandr Sorochynskyi
# On 13 06 18

library(cubature)
library(copula)

source('~/IdV/R/open_neurons.R')

########################################################################################
# utility

ncopula_constrain_param <- function(param) {
  2/(1+exp(-param)) - 1
}

ncopula_release_param <- function(param) {
  -log(2/(param+1) - 1)
}

max_float <- .Machine$double.xmax;

# From copula package source: fitCopula.R
asFinite <- function(x) { # as finite - vectorized in 'x';  NA/NaN basically unchanged
  if(any(nifi <- !is.finite(x)))
    x[nifi] <- sign(x[nifi]) * max_float
  x
}

# the pseudo density for a copula for dim >= 2
generalized_pseudo_density <- function(x, Fu, copula) {
  u <- mapply(function(x,f) f(x), x=x, f=Fu)
  u_minus <- mapply(function(x,f) f(x), x=(x-1), f=Fu)
  ret <- prob(copula,u_minus,u);
  if (ret < 0) warning("Negative density over hyper cube.");
  ret;
  
  # All of this can be done from the library, and faster ;(
  # pd <- function(u) {
  #   pCopula(u, copula);
  # }
  # dim <- length(x);
  # u <- mapply(function(x,f) f(x), x=x, f=Fu)
  # u_minus <- mapply(function(x,f) f(x), x=(x-1), f=Fu)
  # 
  # gen_vec <-  function(index) { tmp <- u; tmp[index] <- u_minus[index]; tmp; }
  # 
  # ret <- sapply(1:dim,
  #              function(i) {
  #                sum(pd(t(combn(dim, i, gen_vec))));
  #              }
  # );
  # # vertecies odd number of steps from origin are to be substracted
  # odd_indecies <- seq(1,length(ret),2);
  # ret[odd_indecies] <- -1*ret[odd_indecies] 
  # # add the origin (dim choose 0)
  # sum(ret) + pd(u)
}

generalized_pseudo_density_ingegration <- function(x, Fu, copula) {
  u <- mapply(function(x,f) f(x), x=x, f=Fu)
  u_minus <- mapply(function(x,f) f(x), x=(x-1), f=Fu)
  
  hcubature(function(t) { dCopula(t,copula) },
            lowerLimit=u_minus,
            upperLimit=u)  
}

make_optim_generalized_pseudo_density <- function(range, dim, cop, Fu) {
  m <- array(NA, dim= rep(range,dim))
  function(a) {
    if (any(a >= range)) {
      return(generalized_pseudo_density(a,Fu,cop));
    }
    index <- matrix(a+1,ncol=dim)
    if (is.na(m[index])) {
      m[index] <<- generalized_pseudo_density(a,Fu,cop)
    }
    m[index]
  }
}

ncopula_make_loglik <- function(x,marginals,copula) {
  function(theta) {
    if (is(copula,"ellipCopula")) {
      theta <- ncopula_constrain_param(theta);
    }
    copula <- setTheta(copula,theta);
    c.density <- make_optim_generalized_pseudo_density(4, ncol(x), copula, marginals);
    likelihoods <- apply(x,1,c.density);
    sum(log(likelihoods));
  }
}

# Can only handle eliptical copulas for now
# x in M(nXdim)
ncopula_fit <- function(x,copula,marginals,inital_val,lower,upper,hessian=FALSE) {
  dim <- ncol(x);
  param_length <- ifelse(is(copula,"rotCopula"),
                            length(copula@copula@param.lowbnd),
                            length(copula@param.lowbnd));
  
  if (is(copula,"indepCopula")) {
    stop("No need to fit indepCopula")
  }
  if (param_length > 1 && !is(copula,"ellipCopula")) {
    stop("ncopula_fit does not support multiple parameter modes for copulas other than elliptical.")
  }
  if (is(copula,"ellipCopula")) {
    optim.method = "Nelder-Mead";
  }
  else {
    optim.method = "Brent";
  }
  
  if (missing(marginals)) marginals <- apply(x,2,ecdf);
  if (missing(inital_val)) {
    if (is(copula,"ellipCopula")) {
      # density becomes too sparse, and becomes zero, numerically
      # indep case seems to be more reliable
      inital_val = rep(0, param_length);
      # if (dim > 4) {
      # }
      # else {
      #   inital_val <- ncopula_release_param(P2p(cov(x)));
      # }
    }
    else {
      inital_val <- runif(param_length)
    }
  }
  if (missing(lower)) {
      if (is(copula, "ellipCopula")) {
          lower = -Inf;
      }
      else {
        if (is(copula,"rotCopula")) {
          lower = rep(max(-20, copula@copula@param.lowbnd),param_length)
        }
        else {
          lower = rep(max(-20, copula@param.lowbnd),param_length)
        }
      }
  }
  if (missing(upper)) {
    if (is(copula, "ellipCopula")) {
      upper = Inf;
    }
    else {
      if (is(copula, "rotCopula")) {
        upper = rep(min(20, copula@copula@param.upbnd),param_length);
      }
      else {
        upper = rep(min(20, copula@param.upbnd),param_length);
      }
    }
  }
  
  loglik <- ncopula_make_loglik(x,marginals,copula);
  
  fit <- optim(inital_val,
               function(x) -loglik(x),
               method = optim.method,
               lower = lower,
               upper=upper,
               hessian= hessian)
  
  if (fit$convergence) {
    warning("Optimization may not have converged: non zero convergencce code from optim()",
              paste("\noptim(): ", fit$message));
  }
  
  if(is(copula,"ellipCopula")) {
    par <- ncopula_constrain_param(fit$par)
  }
  else {
    par <- fit$par
  }
  
  list(par= par,
       hessian = fit$hessian,
       loglik_max= -fit$value); 
}
