# Created by Oleksandr Sorochynskyi
# On 19 07 18

library(bbmle)

source('~/IdV/R/open_neurons.R')
source("copula.R")

################################################################################
# utility

# Generated random data of same dimentions as m
# Using m's empirical marginals
cond_model2.2_rand <- function(n1, n2, m, param, n.rep) {
  ret <- sapply(1:dim(m)[2], function(t) {
    q1 <- make_quantile(m[n1,t,]);
    q2 <- make_quantile(m[n2,t,]);
    rand <- copula_rand_emp(n.rep, frankCopula(param),q1,q2);
    rand;
  },simplify = "array");
  aperm(ret, c(2,3,1))
}

cond_model3_rand <- function(n1, n2, m, param, n.rep) {
  rand <- copula_rand_emp(n.rep*dim(m)[2],frankCopula(param),x1=as.vector(m[n1,,]),x2=as.vector(m[n2,,]));
  ret <- array(NA,dim=c(2, dim(m)[2], n.rep));
  ret[1,,] <- array(rand[,1],dim=dim(ret)[2:3])
  ret[2,,] <- array(rand[,2],dim=dim(ret)[2:3])
  ret;
}

################################################################################
# Fit all 3 options of copulas


# check if we even try to model this bin.
cond_joint_spike_criteria <- function(x,y) {
  any(x*y != 0);
}

cond_mean_prod_criteria <- function(x,y) {
  cond_mean_prod_coef(x,y) != 0;
}

cond_select_bins <- function(n1, n2, m, criteria) {
  ret <- sapply(1:dim(m)[2], function(t) {
    if (criteria(m[n1,t,],m[n2,t,])) {
      t;
    }
    else {
      NA;
    }
  });
  ret[!is.na(ret)];
}

# Option 1: Model only time bins where the neurons spike.
# Fit all bins that satisfy bin_criteria
cond_model1 <- function(n1, n2, m, copula, bin_subset) {
  fit_bin <- function(n1,n2,copula) {
    lapply(bin_subset, function(t) {
      fit_copula(m[n1,t,],m[n2,t,],copula);  
    });
  }
  
  extract_params <- function(fit) {
    params <- sapply(seq_along(fit), function(i) {
      fit[[i]]$par;
    })
  }
  
  extract_loglik <- function(fit) {
    params <- sapply(seq_along(fit), function(i) {
      fit[[i]]$loglik_max;
    })
  }
  
  # Try to fit some cops
  fit <- fit_bin(n1,n2,copula);
  list(params=extract_params(fit),
       loglikes= extract_loglik(fit))
}


# Transform the data
# Estimate the cdf seperatly for every time bin.
# And trasform the data with it.
cond_prob_transform_by_bin <- function(n1, n2, m) {
  ret <- sapply(1:dim(m)[2], function(t) {
    cdf1 <- ecdf(m[n1,t,]); 
    cdf2 <- ecdf(m[n2,t,]); 
    u1 <- cdf1(m[n1,t,]);
    u2 <- cdf2(m[n2,t,]);
    cbind(u1,u2);
  }, simplify="array");
  dim(ret) <- c(dim(m)[3], 2, dim(m)[2])
  ret;
}


# Option 2: Here we suppose that the copula is constant over time, but not the 
#           marginal distributions. So we estimate marginals

# Option 2.1
# Assume that we can trhow all observations into one big pile
cond_model2.1 <- function(n1, n2, m, copula, bin_subset) {
  u <- cond_prob_transform_by_bin(n1,n2,m[,,bin_subset]);
  
  # collapse rep, and time dimentions
  u1 <- as.vector(u[,1,]);
  u2 <- as.vector(u[,2,]);
  
  # do one big fit
  fit <- fit_copula(u1, u2, copula, identity, identity);
  list(param= fit$par,
       loglik= fit$loglik_max)
}

# Option 2.2
# Make a loglikelihood for every time bin, then sum them all up
cond_model2.2 <- function(n1,n2,m,copula,bin_subset) {
  # Transform the data
  # Estimate the cdf seperatly for every time bin.
  # And trasform the data with it.
  
  # So that m doesn't simplyfy to a matrix if dim(bin_subset) == 1
  dim_m <- c(dim(m)[1], dim(bin_subset), dim(m)[3])
  m <- m[,bin_subset,];
  dim(m) <- dim_m
  
  u <- cond_prob_transform_by_bin(n1,n2,m);
  
  # Make a loglikelihood for every time bin, then sum them all up
  # No way to do this using fit_copula
  
  loglikelihood <- function(theta) {
    # Transform theta into parameter bounds:
    #if (is(copula, "normalCopula")) theta <- 2*(1/(exp(-theta)+1))-1 
    
    cop <- setTheta(copula, theta);
    loglike <- sapply(1:dim(u)[3], function(t) {
      pd <- make_pseudo_density(cop, unique(c(0,1,u[,1,t])), unique(c(0,1,u[,2,t])));
      likelihood <- sapply(1:dim(u)[1], function(r) {
        pd(u[r,1,t], u[r,2,t]); 
      });
      sum(log(likelihood));
    });
    -sum(loglike);
  }
  
  fit <- mle2(minuslogl = loglikelihood,
              start= list(theta = 0),
              skip.hessian = TRUE,
              upper=min(20, copula@param.upbnd),
              lower=max(-20, copula@param.lowbnd),
              method="Brent");
  list(param=coef(fit), loglike= logLik(fit))
}


# Option 3: distribution is constant over time bins, and does not depend on S
#       Trhow all obs in one big pile
cond_model3 <- function(n1, n2, m, copula, bin_subset) {
  x <- as.vector(m[n1,,]);
  y <- as.vector(m[n2,,]);
  
  
  fit <- fit_copula(x,y,copula);
  list(param= fit$par,
       loglik= fit$loglik_max)

}

################################################################################
# Mean product coefficient

cond_mean_prod_coef <- function(x,y) {
  sqrt(mean(x)*mean(y));
}

################################################################################
# Covariance decomposition

cond_noise_cov <- function(i,j, m = get_rep_neurons()) {
  mean(
    sapply(1:dim(m)[2], function(t) {
      cov(m[i,t,],m[j,t,]);
    })
  )
}

cond_stimul_cov <- function(i,j, m = get_rep_neurons()) {
  cov(
    sapply(1:dim(m)[2], function(t) {
      mean(m[i,t,])
    })
    ,
    sapply(1:dim(m)[2], function(t) {
      mean(m[j,t,])
    })
  )
}

cond_total_cov <- function(i,j, m = get_rep_neurons()) {
  mean(
    mean(
      sapply(1:dim(m)[2], function(t) {
        mean((m[i,t,] - mean(m[i,,]))*(m[j,t,] - mean(m[j,,])));
      })
    )
  ) 
}
