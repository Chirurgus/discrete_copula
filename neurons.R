# Created by Oleksandr Sorochynskyi
# On 14 06 18

source("~/IdV/R/open_neurons.R");
source("~/IdV/R/copula.R");

neuron.k_fold <- function(x1,x2,copula,k=4,...) {
  res_par <- rep(NA, k);
  res_ll <- rep(NA, k);
  
  index <- sample(seq_along(x1),length(x1),replace=FALSE);
  group <- cut(index, k, labels=FALSE);
  for (i in 1:k) {
    train_i <- index[group == i];
    test_i <- index[group != i]; 
    fit <- fit_copula(x1[train_i],x2[train_i],copula,...);
    test_copula <- setTheta(copula, fit$par);
    ll <- copula_loglik(x1[test_i],x2[test_i],test_copula);
    
    res_par[i] <- fit$par;
    res_ll[i] <- ll
  }
  list(params=res_par, loglik=res_ll);
}

spearman_cor_matrix <- function(m) {
  nrow <- nrow(m)
  cor_m <- matrix(NA, nrow=nrow, ncol=nrow)
  for (i in 1:nrow) {
    for (j in 1:i) {
      cor_m[i,j] <- cor(m[i,],m[j,],method="spearman")
    }
  }
  cor_m;
}

pairwise_copula_fit <- function(m, cop.list, cop.options) {
  stopifnot(length(cop.list) == length(cop.options));
  
  nrow <- nrow(m);
  param <- list()
  
  pb <- txtProgressBar(0, nrow*(nrow-1)/2, style=3)
  
  for (k in 1:length(cop.list)) {
    param[[k]] <- matrix(NA, nrow=nrow, ncol=nrow);
    name(param[[k]]) <- class(cop.list[[k]])
    
    options <- cop.options[[k]];
    options[["copula"]] <- cop.list[[k]];
    for (i in 2:nrow) {
      for (j in 1:(i-1)) {
          options[["x1"]] <- m[i,];
          options[["x2"]] <- m[j,];
          param[[k]][i,j] <- do.call(fit_copula,options)$par;
        }
        
        setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      }
  }

  close(pb);
  
  param;
}

pairwise_loglik <- function(m, cop.list, params) {
  stopifnot(length(cop.list) == length(params));
  
  nrow <- nrow(m);
  loglikes <- rep(NA, length(cop.list));
  
  pb <- txtProgressBar(0, nrow*(nrow-1)/2, style=3)
  
  for (k in 1:length(cop.list)) {
    loglikes[[k]] <- NA;
    name(loglikes[[k]]) <- paste(class(cop.list[[k]]),"_loglike");
    
    for (i in 2:nrow) {
      for (j in 1:(i-1)) {
        loglikes[[k]] <- copula_loglik(m[i,], m[j,]);
        }
        
        setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      }
  }

  close(pb);
  
  param;

}