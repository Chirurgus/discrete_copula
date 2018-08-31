# Created by Oleksandr Sorochynskyi
# On 29 08 18

tri <- function(m) {
  m[lower.tri(m)]
}

pair_near <- c(23,12)
pair_far <- c(23,3)
pair_median <- c(23,2)

# params
load("summary/params.RData")
# covar 
load("summary/covar.RData")
# distances
load("summary/dist.RData")


# 1) infered param: distance param_check param_fullfieldd param_bar
infered_param <- data.frame(distance= tri(distance[1:25,1:25]),
           param_check= tri(params$model2.2$par$check[1:25,1:25]),
           param_fullfield= tri(params$model2.2$par$fullfield[1:25,1:25]),
           param_bar= tri(params$model2.2$par$bar[1:25,1:25]));
write.table(infered_param, "export/time_model/params.txt", row.names=FALSE)

# 2) time-behavior of noise-cov 12&23: noise-cov-data(time), noise-cov-model(time)
data <- bin_checkerboard(pair_near);
rand_check <- cond_model2.2_rand(1, 2, data, params$model2.2$par$check[pair_near[1], pair_near[2]], 1000);

data_noise_cov <- sapply(1:dim(data)[2], function(t) {
      cov(data[1,t,],data[2,t,]);
    });
rand_check_noise_cov <- sapply(1:dim(rand_check)[2], function(t) {
      cov(rand_check[1,t,],rand_check[2,t,]);
    });

time_behavior <- data.frame(noise_check_data_check_param_time_data= data_noise_cov,
                            noise_check_data_check_param_time_model= rand_check_noise_cov);
write.table(time_behavior, "export/time_model/cov_noise_check_data_check_param_time_12_23.txt")

# 3) for (check, fullfield,bars): noise-cov-data noise-cov-model
check <- data.frame(noise_check_data_check_param_data= tri(covar$check_data$check_param$model2.2$real$noise[1:25,1:25]),
                    noise_check_data_check_param_model= tri(covar$check_data$check_param$model2.2$simu$noise[1:25,1:25]));
fullfield <- data.frame(noise_fullfield_data_fullfield_param_data= tri(covar$fullfield_data$fullfield_param$model2.2$real$noise[1:25,1:25]),
                    noise_fullfield_data_fullfield_param_model= tri(covar$fullfield_data$fullfield_param$model2.2$simu$noise[1:25,1:25]));
bar <- data.frame(noise_bar_data_bar_param_data= tri(covar$bar_data$bar_param$model2.2$real$noise[1:25,1:25]),
                    noise_bar_data_bar_param_model= tri(covar$bar_data$bar_param$model2.2$simu$noise[1:25,1:25]));

write.table(check, "export/time_model/cov_noise_check_data_check_param.txt", row.names= FALSE)
write.table(fullfield, "export/time_model/cov_noise_fullfield_data_fullfield_param.txt", row.names= FALSE)
write.table(bar, "export/time_model/cov_noise_bar_data_bar_param.txt", row.names= FALSE)

# 4) with bar_data: noise-cov-data(time), noise-cov-model(time,twobar-param),  noise-cov-model(time,check-param)
data <- bin_barmovie(pair_near);
rand_bar <- cond_model2.2_rand(1, 2, data, params$model2.2$par$bar[pair_near[1], pair_near[2]], 1000);
rand_check <- cond_model2.2_rand(1, 2, data, params$model2.2$par$check[pair_near[1], pair_near[2]], 1000);
rand_fullfield <- cond_model2.2_rand(1, 2, data, params$model2.2$par$fullfield[pair_near[1], pair_near[2]], 1000);

data_noise_cov <- sapply(1:dim(data)[2], function(t) {
      cov(data[1,t,],data[2,t,]);
    });
rand_bar_noise_cov <- sapply(1:dim(rand_bar)[2], function(t) {
      cov(rand_bar[1,t,],rand_bar[2,t,]);
    });
rand_check_noise_cov <- sapply(1:dim(rand_check)[2], function(t) {
      cov(rand_check[1,t,],rand_check[2,t,]);
    });
rand_fullfield_noise_cov <- sapply(1:dim(rand_fullfield)[2], function(t) {
      cov(rand_fullfield[1,t,],rand_fullfield[2,t,]);
    });

time_behavior_bar <- data.frame(noise_bar_data_bar_param_time_data = data_noise_cov,
                                noise_bar_data_bar_param_time_model= rand_bar_noise_cov,
                                noise_bar_data_fullfield_param_time_model= rand_fullfield_noise_cov,
                                noise_bar_data_check_param_time_model= rand_check_noise_cov);
write.table(time_behavior_bar, "export/time_model/cov_noise_bar_data_all_param_time_12_23.txt",row.names= FALSE)

# 5) for (fullfield,bars): noise-cov-data(check_param) noise-cov-model(check_param)
fullfield <- data.frame(noise_fullfield_data_fullfield_param_data= tri(covar$fullfield_data$check_param$model2.2$real$noise[1:25,1:25]),
                    noise_fullfield_data_fullfield_param_model= tri(covar$fullfield_data$check_param$model2.2$simu$noise[1:25,1:25]));
bar <- data.frame(noise_bar_data_bar_param_data= tri(covar$bar_data$check_param$model2.2$real$noise[1:25,1:25]),
                    noise_bar_data_bar_param_model= tri(covar$bar_data$check_param$model2.2$simu$noise[1:25,1:25]));

write.table(fullfield, "export/time_model/cov_noise_fullfield_data_check_param.txt", row.names= FALSE)
write.table(bar, "export/time_model/cov_noise_bar_data_check_param.txt", row.names= FALSE)