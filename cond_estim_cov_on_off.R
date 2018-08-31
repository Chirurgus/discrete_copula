# Created by Oleksandr Sorochynskyi
# On 17 08 18

source('~/IdV/R/conditional.R')

n0 <- 25
n1 <- 10
m0 <- bin_fullfield(1:n0,cell.type= 1);
m1 <- bin_fullfield(1:n1,cell.type= 6);
m <- array(NA, dim=c(n1+n0,dim(m0)[2],dim(m0)[3]))
m[1:n0,,] <- m0;
m[-(1:n0),,] <- m1;
n.rep <- 1000;
filename <- "cov/conditional_all_on_off_bar_cov_bar_param_1000rep.RData";

param_filename <- "params/conditional_all_on_off_bar_frank_nonzero_mean.RData";
load(param_filename);

noise_cov_real <- matrix(NA, ncol=n1, nrow=n0);
noise_cov_simu <- matrix(NA, ncol=n1, nrow=n0);
stimul_cov_real <- matrix(NA, ncol=n1, nrow=n0);
stimul_cov_simu <- matrix(NA, ncol=n1, nrow=n0);

pb <- txtProgressBar(min=0, n0*n1*2, style= 3);

for (i in 1:n0) {
  for (j in 1:n1) {
    rand <- cond_model2.2_rand(i, n0+j, m, model2.2$param[i,j], n.rep);
    
    noise_cov_real[i,j] <- cond_noise_cov(i,j+n0, m);
    stimul_cov_real[i,j] <- cond_stimul_cov(i, n0+j, m);
    
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
noise_cov_real <- matrix(NA, ncol=n1, nrow=n0);
noise_cov_simu <- matrix(NA, ncol=n1, nrow=n0);
stimul_cov_real <- matrix(NA, ncol=n1, nrow=n0);
stimul_cov_simu <- matrix(NA, ncol=n1, nrow=n0);



for (i in 1:n0) {
  for (j in 1:n1) {
    rand <- cond_model3_rand(i, n0+j, m, model3$param[i,j], n.rep);
    
    noise_cov_real[i,j] <- cond_noise_cov(i,j+n0, m);
    stimul_cov_real[i,j] <- cond_stimul_cov(i, n0+j, m);
    
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
