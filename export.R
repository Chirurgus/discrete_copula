# Created by Oleksandr Sorochynskyi 
# On 26 08 18

library(R.matlab)

export_params <- function(file.in, file.out) {
  load(file.in);
  writeMat(con=file.out,
           time_model_param= model2.2$param,
           time_model_loglik= model2.2$loglik)  
}

# All paramete file names
files.in <- c("params/conditional_1_25_off_bar_frank_nonzero_mean.RData",
              "params/conditional_1_25_off_check_frank_nonzero_mean.RData",
              "params/conditional_1_25_off_fullfield_frank_nonzero_mean.RData",
              "params/conditional_1_19_on_bar_frank_nonzero_mean.RData",
              "params/conditional_1_19_on_check_frank_nonzero_mean.RData",
              "params/conditional_1_19_on_fullfield_frank_nonzero_mean.RData",
              "params/conditional_all_on_off_bar_frank_nonzero_mean.RData",
              "params/conditional_all_on_off_check_frank_nonzero_mean.RData",
              "params/conditional_all_on_off_fullfield_frank_nonzero_mean.RData");

files.out <- c("export/off_bar.mat",
               "export/off_check.mat",
               "export/off_fullfield.mat",
               "export/on_bar.mat",
               "export/on_check.mat",
               "export/on_fullfield.mat",
               "export/on_off_bar.mat",
               "export/on_off_check.mat",
               "export/on_off_fullfield.mat");

mapply(export_params, files.in, files.out)






#######
pair <- c(12,23)
n.rep <- 1000;

load("params/conditional_1_25_off_check_frank_nonzero_mean.RData");
check_param <- model2.2$param[pair[2],pair[1]];
load("params/conditional_1_25_off_bar_frank_nonzero_mean.RData");
bar_param <- model2.2$param[pair[2],pair[1]];
load("params/conditional_1_25_off_fullfield_frank_nonzero_mean.RData");
fullfield_param <- model2.2$param[pair[2],pair[1]];

data <- bin_checkerboard(pair);
rand_check <- cond_model2.2_rand(1, 2, data, check_param, n.rep);
#rand_bar <- cond_model2.2_rand(1, 2, data, bar_param, n.rep);
#rand_fullfield <- cond_model2.2_rand(1, 2, data, fullfield_param, n.rep);

data_noise_cov <- sapply(1:dim(data)[2], function(t) {
      cov(data[1,t,],data[2,t,]);
    });
rand_check_noise_cov <- sapply(1:dim(rand_check)[2], function(t) {
      cov(rand_check[1,t,],rand_check[2,t,]);
    });
#rand_bar_noise_cov <- sapply(1:dim(rand_bar)[2], function(t) { cov(rand_bar[1,t,],rand_bar[2,t,]); })
#rand_fullfield_noise_cov <- sapply(1:dim(rand_fullfield)[2], function(t) { cov(rand_fullfield[1,t,],rand_fullfield[2,t,]); })

writeMat(con="export/12_23_off_check_noise_cov_check_data.mat",
         check_data_noise_cov=data_noise_cov,
         simu_check_noise_cov=rand_check_noise_cov);
         #simu_fullfield_noise_cov=rand_fullfield_noise_cov,
         #simu_bar_noise_cov=rand_bar_noise_cov);


#####

export_cov <- function(file.in, file.out) {
  load(file.in);
  writeMat(con=file.out,
           data_noise_cov= model2.2_cov_real$noise,
           simu_noise_cov= model2.2_cov_simu$noise);
           #data_stimul_cov= model2.2_cov_real$stimul);
           #simu_stimul_cov= model2.2_cov_simu$stimul);
}

files.in <- c("cov/conditional_1_25_off_check_cov_check_param_1000rep.RData",
              "cov/conditional_1_25_off_bar_cov_bar_param_1000rep.RData",
              "cov/conditional_1_25_off_fullfield_cov_fullfield_param_1000rep.RData")

files.out <- c("export/off_check_cov_check_param.mat",
               "export/off_bar_cov_bar_param.mat",
               "export/off_fullfield_cov_fullfield_param.mat")

mapply(export_cov, files.in, files.out);


#####

pair <- c(12,23)
n.rep <- 1000;

load("params/conditional_1_25_off_check_frank_nonzero_mean.RData");
check_param <- model2.2$param[pair[2],pair[1]];
load("params/conditional_1_25_off_bar_frank_nonzero_mean.RData");
bar_param <- model2.2$param[pair[2],pair[1]];
#load("params/conditional_1_25_off_fullfield_frank_nonzero_mean.RData");
#fullfield_param <- model2.2$param[pair[2],pair[1]];

data <- bin_barmovie(pair);
rand_check <- cond_model2.2_rand(1, 2, data, check_param, n.rep);
rand_bar <- cond_model2.2_rand(1, 2, data, bar_param, n.rep);
#rand_fullfield <- cond_model2.2_rand(1, 2, data, fullfield_param, n.rep);

data_noise_cov <- sapply(1:dim(data)[2], function(t) { cov(data[1,t,],data[2,t,]); })
rand_check_noise_cov <- sapply(1:dim(rand_check)[2], function(t) { cov(rand_check[1,t,],rand_check[2,t,]); })
rand_bar_noise_cov <- sapply(1:dim(rand_bar)[2], function(t) { cov(rand_bar[1,t,],rand_bar[2,t,]); })
#rand_fullfield_noise_cov <- sapply(1:dim(rand_fullfield)[2], function(t) { cov(rand_fullfield[1,t,],rand_fullfield[2,t,]); })

writeMat(con="export/12_23_off_bar_noise_cov_check_bar_param.mat",
         check_data_noise_cov=data_noise_cov,
         simu_check_noise_cov=rand_check_noise_cov,
         #simu_fullfield_noise_cov=rand_fullfield_noise_cov,
         simu_bar_noise_cov=rand_bar_noise_cov);



#######

files.in <- c("cov/conditional_1_25_off_check_cov_check_param_1000rep.RData",
              "cov/conditional_1_25_off_bar_cov_check_param_1000rep.RData",
              "cov/conditional_1_25_off_fullfield_cov_check_param_1000rep.RData")

files.out <- c("export/off_check_cov_check_param.mat",
               "export/off_bar_cov_check_param.mat",
               "export/off_fullfield_cov_check_param.mat")

mapply(export_cov, files.in, files.out);









# #####
# # Synthetic approach (pull all neuron types together)
# 
# concat_matrix <- function(top.left, bottom.left, bottom.right) {
#   stopifnot(ncol(top.left) == ncol(bottom.left),
#             nrow(bottom.left) == nrow(bottom.right));
#   ret <- matrix(NA,
#                 nrow = nrow(top.left) + nrow(bottom.left),
#                 ncol = ncol(bottom.left) + ncol(bottom.right));
#   ret[1:nrow(top.left), 1:ncol(top.left)] <- top.left;
#   ret[-(1:nrow(top.left)), 1:ncol(top.left)] <- bottom.left;
#   ret[-(1:nrow(top.left)),-(1:ncol(top.left))] <- bottom.right;
#   ret;                
# }
# 
# # params
# 
# load("params/conditional_all_on_off_bar_frank_nonzero_mean.RData")
# on_off_3 <- t(model3$param)
# on_off_3_ll <- t(model3$loglik)
# on_off_2.2 <- t(model2.2$param)
# on_off_2.2_ll <- t(model2.2$loglik)
# load("params/conditional_1_25_off_bar_frank_nonzero_mean.RData");
# off_3 <- model3$param
# off_3_ll <- model3$loglik
# off_2.2 <- model2.2$param
# off_2.2_ll <- model2.2$loglik
# load("params/conditional_1_19_on_bar_frank_nonzero_mean.RData")
# on_3 <- model3$param
# on_3_ll <- model3$loglik
# on_2.2 <- model2.2$param
# on_2.2_ll <- model2.2$loglik
# 
# bar_3 <- concat_matrix(off_3, on_off_3, on_3);
# bar_ll_3 <- concat_matrix(off_3_ll, on_off_3_ll, on_3_ll);
# bar_2.2 <- concat_matrix(off_2.2, on_off_2.2, on_2.2);
# bar_ll_2.2 <- concat_matrix(off_2.2_ll, on_off_2.2_ll, on_2.2_ll);
# 
# load("params/conditional_all_on_off_fullfield_frank_nonzero_mean.RData")
# on_off_3 <- t(model3$param)
# on_off_3_ll <- t(model3$loglik)
# on_off_2.2 <- t(model2.2$param)
# on_off_2.2_ll <- t(model2.2$loglik)
# load("params/conditional_1_25_off_fullfield_frank_nonzero_mean.RData");
# off_3 <- model3$param
# off_3_ll <- model3$loglik
# off_2.2 <- model2.2$param
# off_2.2_ll <- model2.2$loglik
# load("params/conditional_1_19_on_fullfield_frank_nonzero_mean.RData")
# on_3 <- model3$param
# on_3_ll <- model3$loglik
# on_2.2 <- model2.2$param
# on_2.2_ll <- model2.2$loglik
# 
# fullfield_3 <- concat_matrix(off_3, on_off_3, on_3);
# fullfield_ll_3 <- concat_matrix(off_3_ll, on_off_3_ll, on_3_ll);
# fullfield_2.2 <- concat_matrix(off_2.2, on_off_2.2, on_2.2);
# fullfield_ll_2.2 <- concat_matrix(off_2.2_ll, on_off_2.2_ll, on_2.2_ll);
# 
# load("params/conditional_all_on_off_check_frank_nonzero_mean.RData")
# on_off_3 <- t(model3$param)
# on_off_3_ll <- t(model3$loglik)
# on_off_2.2 <- t(model2.2$param)
# on_off_2.2_ll <- t(model2.2$loglik)
# load("params/conditional_1_25_off_check_frank_nonzero_mean.RData");
# off_3 <- model3$param
# off_3_ll <- model3$loglik
# off_2.2 <- model2.2$param
# off_2.2_ll <- model2.2$loglik
# load("params/conditional_1_19_on_check_frank_nonzero_mean.RData")
# on_3 <- model3$param
# on_3_ll <- model3$loglik
# on_2.2 <- model2.2$param
# on_2.2_ll <- model2.2$loglik
# 
# check_3 <- concat_matrix(off_3, on_off_3, on_3);
# check_ll_3 <- concat_matrix(off_3_ll, on_off_3_ll, on_3_ll);
# check_2.2 <- concat_matrix(off_2.2, on_off_2.2, on_2.2);
# check_ll_2.2 <- concat_matrix(off_2.2_ll, on_off_2.2_ll, on_2.2_ll);
# 
# 
# 
# model3 <- list(check=check_3, bar=bar_3, fullfield=fullfield_3)
# model3_ll <- list(check=check_ll_3, bar=bar_ll_3, fullfield=fullfield_ll_3)
# model2.2 <- list(check=check_2.2, bar=bar_2.2, fullfield=fullfield_2.2)
# model2.2_ll <- list(check=check_ll_2.2, bar=bar_ll_2.2, fullfield=fullfield_ll_2.2)
# 
# params <- list(model3=list(par=model3,loglik=model3_ll),
#                model2.2=list(par=model2.2,loglik=model2.2_ll));
# save(file="par.RData", par=params);
# 
# 
# # Now distance
# distance <- concat_matrix(distance_off, distance_on_off, distance_on)
# save(file="dist.RData", dist=distance);
# 
# # Now covariance
# 
