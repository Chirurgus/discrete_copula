# Created by Oleksandr Sorochynskyi
# On 27 7 18

# Insert your wd here
#setwd("/mnt/c/Users/Alexander/Desktop");

#install.packages(c("bbmle", "R.matlab", "copula"));
library(bbmle)
library(R.matlab)
library(copula)

source('~/IdV/R/conditional.R')

#system.time({
################################################################################
# Now do the actual estimation
copula <- frankCopula();
#n <- 2
n0 <- 25
n1 <- 19
#n <- 17
#n <- dim(m)[1]
m0 <- bin_fullfield(1:n0,cell.type= 1);
m1 <- bin_fullfield(1:n1,cell.type= 6);
m <- array(NA, dim=c(n1+n0,dim(m0)[2],dim(m0)[3]))
m[1:n0,,] <- m0;
m[-(1:n0),,] <- m1;
#criteria <- cond_joint_spike_criteria
criteria <- cond_mean_prod_criteria
filename <- "params/conditional_all_on_off_fullfield_frank_nonzero_mean.RData";
existing_filename <- "params/conditional_all_on_off_fullfield_frank_nonzero_mean.RData";

# model1_param <- array(NA, dim=c(n,n,dim(m)[2]));
# model2.1_param <- matrix(NA, ncol=n,nrow=n);
model2.2_param <- matrix(NA, ncol=n1,nrow=n0);
model3_param <- matrix(NA, ncol=n1,nrow=n0);

# model1_loglik <- array(NA, dim=c(n,n,dim(m)[2]));
# model2.1_loglik <- matrix(NA, ncol=n,nrow=n);
model2.2_loglik <- matrix(NA, ncol=n1,nrow=n0);
model3_loglik <- matrix(NA, ncol=n1,nrow=n0);

old_n1 <- 10
old_n0 <- 25
load(existing_filename)
model2.2_param[1:old_n0,1:old_n1] <- model2.2$param
model2.2_loglik[1:old_n0,1:old_n1] <- model2.2$loglik
model3_param[1:old_n0,1:old_n1] <- model3$param
model3_loglik[1:old_n0,1:old_n1] <- model3$loglik

pb <- txtProgressBar(min=0, n1*n0, style= 3);

for (i in 1:n0) {
  for (j in 1:n1) {
    if (i <= old_n0 && j <= old_n1) next;
    
    selected.bins <- cond_select_bins(i,n0+ j, m, criteria);
    
    # fit1 <- tryCatch({
    #   cond_model1(i,j,m,copula,selected.bins);
    # },
    # error= function(e) {
    #   warning(paste("Could't fit model 1 for pair:",i,",",j,"."))
    #   list(params= rep(NA, length(selected.bins)),
    #        loglikes= rep(NA, length(selected.bins)));
    # });
    # 
    # fit2.1 <- tryCatch({
    #   cond_model2.1(i,j,m,copula,selected.bins);
    # },
    # error= function(e) {
    #   warning(paste("Could't fit model 2.1 for pair:",i,",",j,"."))
    #   list(param=NA, loglik= NA);
    # });
    # 
    fit2.2 <- tryCatch({
      cond_model2.2(i,n0+j,m,copula,selected.bins);
    },
    error= function(e) {
      warning(paste("Could't fit model 2.2 for pair:",i,",",j,"."))
      list(param=NA, loglik= NA);
    });
    
    fit3 <- tryCatch({
      cond_model3(i,n0+j,m,copula,selected.bins);
    },
    error= function(e) {
      warning(paste("Could't fit model 3 for pair:",i,",",j,"."))
      list(param=NA, loglik= NA);
    });
    
    # model1_param[i,j,selected.bins] <- fit1$params;
    # model1_loglik[i,j,selected.bins] <- fit1$loglikes;
    # model2.1_param[i,j] <- fit2.1$param;
    # model2.1_loglik[i,j] <- fit2.1$loglik;
    model2.2_param[i,j] <- fit2.2$param;
    model2.2_loglik[i,j] <- fit2.2$loglik;
    model3_param[i,j] <- fit3$param;
    model3_loglik[i,j] <- fit3$loglik;
    
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1);
  }
}

# model1 <- list(param= model1_param,
#                loglik= model1_loglik);
# model2.1 <- list(param= model2.1_param,
#                loglik= model2.1_loglik);
model2.2 <- list(param= model2.2_param,
               loglik= model2.2_loglik);
model3 <- list(param= model3_param,
               loglik= model3_loglik);

#}); # system.time

#save(model1, model2.1, model2.2, model3, file=filename);
save(model2.2, model3, file=filename);
#load(filename);

close(pb);




