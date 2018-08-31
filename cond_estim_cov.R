# Created by Oleksandr Sorochynskyi
# On 17 08 18

source('~/IdV/R/conditional.R')

n <- 25
m <- bin_barmovie(1:n,cell.type = 1);
n.rep <- 1000;
filename <- "cov/conditional_1_25_off_bar_cov_check_param_1000rep.RData";

param_filename <- "params/conditional_1_25_off_check_frank_nonzero_mean.RData";
load(param_filename);

noise_cov_real <- matrix(NA, ncol=n, nrow=n);
noise_cov_simu <- matrix(NA, ncol=n, nrow=n);
stimul_cov_real <- matrix(NA, ncol=n, nrow=n);
stimul_cov_simu <- matrix(NA, ncol=n, nrow=n);

pb <- txtProgressBar(min=0, n*(n-1), style= 3);

for (i in 2:n) {
  for (j in 1:(i-1)) {
    rand <- cond_model2.2_rand(i, j, m, model2.2$param[i,j], n.rep);
    
    noise_cov_real[i,j] <- cond_noise_cov(i,j, m);
    stimul_cov_real[i,j] <- cond_stimul_cov(i, j, m);
    
    noise_cov_simu[i,j] <- cond_noise_cov(1,2,rand);
    stimul_cov_simu[i,j] <- cond_stimul_cov(1,2,rand);
    
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1);
  }
}

total_cov_real <- stimul_cov_real + noise_cov_real;
total_cov_simu <- stimul_cov_simu + noise_cov_simu;

model2.2_cov_real <- list(total=total_cov_real,
                        noise=noise_cov_real,
                        stimul=stimul_cov_real)
model2.2_cov_simu <- list(total=total_cov_simu,
                        noise=noise_cov_simu,
                        stimul=stimul_cov_simu)

# Now model3
noise_cov_real <- matrix(NA, ncol=n, nrow=n);
noise_cov_simu <- matrix(NA, ncol=n, nrow=n);
stimul_cov_real <- matrix(NA, ncol=n, nrow=n);
stimul_cov_simu <- matrix(NA, ncol=n, nrow=n);
n

for (i in 2:n) {
  for (j in 1:(i-1)) {
    rand <- cond_model3_rand(i, j, m, model3$param[i,j], n.rep);
    
    noise_cov_real[i,j] <- cond_noise_cov(i,j,m);
    stimul_cov_real[i,j] <- cond_stimul_cov(i, j, m);
    
    noise_cov_simu[i,j] <- cond_noise_cov(1,2,rand);
    stimul_cov_simu[i,j] <- cond_stimul_cov(1,2,rand);
    
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1);
  }
}

close(pb);

total_cov_real <- stimul_cov_real + noise_cov_real;
total_cov_simu <- stimul_cov_simu + noise_cov_simu;

model3_cov_real <- list(total=total_cov_real,
                        noise=noise_cov_real,
                        stimul=stimul_cov_real)
model3_cov_simu <- list(total=total_cov_simu,
                        noise=noise_cov_simu,
                        stimul=stimul_cov_simu)

save(model2.2_cov_real, model2.2_cov_simu, model3_cov_real, model3_cov_simu, file=filename);
